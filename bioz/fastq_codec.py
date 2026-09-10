"""
Format-aware FASTQ codec.

Instead of compressing the interleaved 4-line-per-read text directly (what
gzip/zstd/xz see when pointed at a raw .fastq), this splits a FASTQ file
into four homogeneous streams that each compress far better on their own:

  ids     -- read headers (share a lot of structure/prefix across reads)
  plus    -- the '+' separator lines (near-constant -> compresses to ~nothing)
  seq     -- all bases concatenated, 2-bit packed (A/C/G/T -> 2 bits each,
             any other symbol like N recorded as a small exception list so
             this step is always exactly lossless)
  qual    -- all quality scores concatenated; optionally binned (lossy) to
             shrink the entropy that dominates most fastq.gz files

Each stream is then run through the backend (best-of zstd/xz). Sequence
and read structure (ids, lengths, plus-lines) are ALWAYS lossless,
regardless of the quality-lossy setting -- only quality scores are ever
approximated, and only if the caller opts into it.

Both directions process the file in bounded-size chunks of records and
stream every intermediate stream through temp files rather than Python
lists/bytes objects, so peak memory stays roughly constant regardless of
input file size (the 2-bit packer keeps only a 0-3-base carry between
chunks so the packed output is byte-identical to packing the whole
concatenated sequence at once).
"""

import gzip
import struct
import tempfile
from pathlib import Path

import numpy as np

from . import backend
from . import container

BASE_TO_CODE = {65: 0, 67: 1, 71: 2, 84: 3}  # A C G T (ASCII)
CODE_TO_BASE = bytes([65, 67, 71, 84])

# Vectorized lookup table (numpy), built once: base byte value -> 2-bit code,
# or 255 as a sentinel for "not a plain ACGT base" (N, lowercase, IUPAC, ...).
_BASE_LUT = np.full(256, 255, dtype=np.uint8)
for _b, _c in BASE_TO_CODE.items():
    _BASE_LUT[_b] = _c
_CODE_LUT = np.frombuffer(CODE_TO_BASE, dtype=np.uint8)

N_QUALITY_BINS = 8
QUALITY_BIN_WIDTH = 8  # covers Phred 0..63 in 8 bins of width 8

_EXC_DTYPE = np.dtype([("pos", "<u8"), ("ch", "u1")])  # packed, 9 bytes/exception

CHUNK_READS = 200_000  # records buffered in memory at a time, both directions


def _opener(path):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")


def _parse_fastq(path):
    """Yields (id_line, seq_line, plus_line, qual_line) as bytes, newline-stripped."""
    with _opener(path) as f:
        while True:
            id_line = f.readline()
            if not id_line:
                return
            seq_line = f.readline()
            plus_line = f.readline()
            qual_line = f.readline()
            if not qual_line:
                raise ValueError("truncated FASTQ record (file doesn't end on a 4-line boundary)")
            id_line = id_line.rstrip(b"\n")
            plus_line = plus_line.rstrip(b"\n")
            if not id_line.startswith(b"@"):
                raise ValueError(f"expected '@' header line, got {id_line!r}")
            if not plus_line.startswith(b"+"):
                raise ValueError(f"expected '+' separator line, got {plus_line!r}")
            yield (
                id_line[1:],
                seq_line.rstrip(b"\n"),
                plus_line[1:],
                qual_line.rstrip(b"\n"),
            )


def _bin_quality(qual_concat: bytes) -> bytes:
    # Phred+33 ASCII -> bin index -> representative Phred+33 ASCII, so the
    # output stream is still printable/ASCII-shaped (helps generic backends)
    # but only N_QUALITY_BINS distinct byte values appear, collapsing entropy.
    table = bytes(
        min((q // QUALITY_BIN_WIDTH) * QUALITY_BIN_WIDTH + QUALITY_BIN_WIDTH // 2, 63) + 33
        for q in range(256 - 33)
    )
    # bytes below 33 shouldn't occur in valid Phred+33 quality; map identity as a safety net
    full_table = bytes(range(33)) + table
    return qual_concat.translate(full_table)


def _translate_and_log_exceptions(chunk_arr, base_offset, exc_f):
    codes = _BASE_LUT[chunk_arr].copy()
    exc_mask = codes == 255
    exc_positions = np.nonzero(exc_mask)[0]
    if exc_positions.size:
        exc_arr = np.zeros(exc_positions.size, dtype=_EXC_DTYPE)
        exc_arr["pos"] = (exc_positions + base_offset).astype("<u8")
        exc_arr["ch"] = chunk_arr[exc_positions]
        exc_f.write(exc_arr.tobytes())
        codes[exc_mask] = 0
    return codes


def _pack_codes(codes, out_f):
    """codes length must be a multiple of 4."""
    codes = codes.reshape(-1, 4)
    packed = (codes[:, 0] | (codes[:, 1] << 2) | (codes[:, 2] << 4) | (codes[:, 3] << 6)).astype(np.uint8)
    out_f.write(packed.tobytes())


class _StreamingSeqPacker:
    """2-bit packs sequence data fed in arbitrary-sized pieces, producing a
    byte stream identical to packing the whole concatenation at once --
    keeps only a 0-3-base carry (translated codes, not raw bytes) between
    calls, so exceptions are logged from real bases only, never padding."""

    def __init__(self, out_f, exc_f):
        self.out_f = out_f
        self.exc_f = exc_f
        self.carry_codes = np.zeros(0, dtype=np.uint8)
        self.global_pos = 0

    def feed(self, seq_bytes: bytes):
        if not seq_bytes:
            return
        arr = np.frombuffer(seq_bytes, dtype=np.uint8)
        codes = _translate_and_log_exceptions(arr, self.global_pos, self.exc_f)
        self.global_pos += arr.size
        combined = np.concatenate([self.carry_codes, codes]) if self.carry_codes.size else codes
        n_full = (combined.size // 4) * 4
        if n_full:
            _pack_codes(combined[:n_full], self.out_f)
        self.carry_codes = combined[n_full:]

    def finish(self):
        if self.carry_codes.size:
            pad = (-self.carry_codes.size) % 4
            padded = np.concatenate([self.carry_codes, np.zeros(pad, dtype=np.uint8)]) if pad else self.carry_codes
            _pack_codes(padded, self.out_f)
            self.carry_codes = np.zeros(0, dtype=np.uint8)


def compress_file(in_path, out_path, threads: int = 0, try_xz: bool = True, lossy_quality: bool = False, ultra: bool = False):
    tmp_dir = tempfile.mkdtemp(prefix="bioz_fastq_")
    tmp_dir = Path(tmp_dir)
    paths = {name: tmp_dir / name for name in ("ids", "plus", "seq_packed", "seq_exc", "lengths", "qual")}
    orig_size = 0
    n_reads = 0

    try:
        with open(paths["ids"], "wb") as ids_f, \
             open(paths["plus"], "wb") as plus_f, \
             open(paths["seq_packed"], "wb") as seqp_f, \
             open(paths["seq_exc"], "wb") as exc_f, \
             open(paths["lengths"], "wb") as len_f, \
             open(paths["qual"], "wb") as qual_f:

            packer = _StreamingSeqPacker(seqp_f, exc_f)
            qual_buf = []
            qual_buf_size = 0

            for id_line, seq_line, plus_line, qual_line in _parse_fastq(in_path):
                if len(seq_line) != len(qual_line):
                    raise ValueError(f"sequence/quality length mismatch for read {id_line!r}")
                ids_f.write(id_line + b"\n")
                plus_f.write(plus_line + b"\n")
                len_f.write(struct.pack("<I", len(seq_line)))
                packer.feed(seq_line)
                qual_buf.append(qual_line)
                qual_buf_size += len(qual_line)
                orig_size += len(id_line) + 1 + len(seq_line) + 1 + len(plus_line) + 1 + len(qual_line) + 1
                n_reads += 1

                if qual_buf_size >= (1 << 20):
                    chunk = b"".join(qual_buf)
                    qual_f.write(_bin_quality(chunk) if lossy_quality else chunk)
                    qual_buf, qual_buf_size = [], 0

            if qual_buf:
                chunk = b"".join(qual_buf)
                qual_f.write(_bin_quality(chunk) if lossy_quality else chunk)
            packer.finish()

        total_bases = packer.global_pos

        streams_meta = [(paths["ids"], False), (paths["plus"], False), (paths["seq_packed"], False),
                         (paths["seq_exc"], False), (paths["lengths"], False), (paths["qual"], False)]
        streams = []
        for path, _ in streams_meta:
            backend_id, comp_path, _size = backend.compress_file(path, threads=threads, try_xz=try_xz, ultra=ultra)
            streams.append((backend_id, comp_path))

        meta = struct.pack("<BQQ", 1 if lossy_quality else 0, n_reads, total_bases)
        meta_backend_id, meta_blob = backend.compress_bytes(meta, threads=threads, try_xz=False)

        container.write_container(
            out_path, container.FORMAT_FASTQ, [(meta_backend_id, meta_blob)] + streams
        )
        for backend_id, comp_path in streams:
            Path(comp_path).unlink(missing_ok=True)

        comp_size = Path(out_path).stat().st_size
        return orig_size, comp_size
    finally:
        for p in paths.values():
            p.unlink(missing_ok=True)
        tmp_dir.rmdir()


class _StreamingSeqUnpacker:
    """Reverses _StreamingSeqPacker: hands out `n` bases at a time (matching
    each read's stored length), reading the packed stream and applying
    exceptions lazily in forward order -- never materializes the whole
    decoded sequence at once."""

    def __init__(self, seqp_path, exc_path, total_bases, read_chunk_bases=4 << 20):
        self._seqp_f = open(seqp_path, "rb")
        self._exc_f = open(exc_path, "rb")
        self._total_bases = total_bases
        self._read_chunk_bases = read_chunk_bases  # must be a multiple of 4
        self._buf = np.zeros(0, dtype=np.uint8)  # decoded, not-yet-consumed bases
        self._buf_start_pos = 0  # absolute position of self._buf[0]
        self._next_exc = self._exc_f.read(9)

    def _refill(self):
        packed = self._seqp_f.read(self._read_chunk_bases // 4)
        if not packed:
            return False
        packed_arr = np.frombuffer(packed, dtype=np.uint8)
        codes = np.empty(packed_arr.size * 4, dtype=np.uint8)
        codes[0::4] = packed_arr & 0x3
        codes[1::4] = (packed_arr >> 2) & 0x3
        codes[2::4] = (packed_arr >> 4) & 0x3
        codes[3::4] = (packed_arr >> 6) & 0x3
        bases = _CODE_LUT[codes].copy()

        start = self._buf_start_pos + self._buf.size
        while self._next_exc:
            (pos,) = struct.unpack("<Q", self._next_exc[:8])
            if pos >= start + bases.size:
                break
            ch = self._next_exc[8]
            bases[pos - start] = ch
            self._next_exc = self._exc_f.read(9)

        self._buf = np.concatenate([self._buf, bases]) if self._buf.size else bases
        return True

    def read_n(self, n: int) -> bytes:
        while self._buf.size < n:
            if not self._refill():
                break
        out = self._buf[:n].tobytes()
        self._buf = self._buf[n:]
        self._buf_start_pos += n
        return out

    def close(self):
        self._seqp_f.close()
        self._exc_f.close()


def decompress_file(in_path, out_path):
    format_id, index = container.iter_container_streams(in_path)
    if format_id != container.FORMAT_FASTQ:
        raise ValueError("not a fastq-format .bioz file")

    (meta_backend, meta_off, meta_len), *rest_idx = index
    tmp_dir = Path(tempfile.mkdtemp(prefix="bioz_fastq_dec_"))
    meta_c = tmp_dir / "meta.c"
    try:
        container.extract_stream_to_file(in_path, meta_off, meta_len, meta_c)
        meta = backend.decompress_bytes(*_read_small(meta_c, meta_backend))
        lossy_quality, n_reads, total_bases = struct.unpack("<BQQ", meta)

        names = ["ids", "plus", "seq_packed", "seq_exc", "lengths", "qual"]
        extracted = {}
        decompressed = {}
        for name, (b_id, off, length) in zip(names, rest_idx):
            ex_path = tmp_dir / f"{name}.c"
            container.extract_stream_to_file(in_path, off, length, ex_path)
            dec_path = tmp_dir / name
            backend.decompress_to_file(b_id, ex_path, dec_path)
            ex_path.unlink(missing_ok=True)
            decompressed[name] = dec_path

        ids_f = open(decompressed["ids"], "rb")
        plus_f = open(decompressed["plus"], "rb")
        len_f = open(decompressed["lengths"], "rb")
        qual_f = open(decompressed["qual"], "rb")
        seq_unpacker = _StreamingSeqUnpacker(decompressed["seq_packed"], decompressed["seq_exc"], total_bases)

        try:
            with open(out_path, "wb") as out_f:
                for _ in range(n_reads):
                    id_line = ids_f.readline().rstrip(b"\n")
                    plus_line = plus_f.readline().rstrip(b"\n")
                    (length,) = struct.unpack("<I", len_f.read(4))
                    seq = seq_unpacker.read_n(length)
                    qual = qual_f.read(length)
                    out_f.write(b"@" + id_line + b"\n")
                    out_f.write(seq + b"\n")
                    out_f.write(b"+" + plus_line + b"\n")
                    out_f.write(qual + b"\n")
        finally:
            ids_f.close()
            plus_f.close()
            len_f.close()
            qual_f.close()
            seq_unpacker.close()

        return bool(lossy_quality)
    finally:
        for p in tmp_dir.iterdir():
            p.unlink(missing_ok=True)
        tmp_dir.rmdir()


def _read_small(path, backend_id):
    """Helper: (backend_id, bytes) for a small already-extracted stream file,
    matching backend.decompress_bytes's signature."""
    return backend_id, Path(path).read_bytes()
