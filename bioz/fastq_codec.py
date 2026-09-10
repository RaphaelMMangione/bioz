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

Each stream is then run through backend.compress_bytes (best-of zstd/xz).
Sequence and read structure (ids, lengths, plus-lines) are ALWAYS lossless,
regardless of the quality-lossy setting -- only quality scores are ever
approximated, and only if the caller opts into it.
"""

import gzip
import struct

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


_EXC_DTYPE = np.dtype([("pos", "<u8"), ("ch", "u1")])  # packed, 9 bytes/exception


def _pack_seq(seq_concat: bytes):
    """2-bit pack ACGT (vectorized with numpy); return (packed_bytes, exceptions_bytes).
    Any byte that isn't a plain uppercase A/C/G/T (N, lowercase, IUPAC codes, ...)
    is recorded as a (position, original_byte) exception rather than lossily
    coerced, so this step is exactly lossless regardless of input alphabet."""
    n = len(seq_concat)
    arr = np.frombuffer(seq_concat, dtype=np.uint8)
    codes = _BASE_LUT[arr].copy()

    exc_mask = codes == 255
    exc_positions = np.nonzero(exc_mask)[0]
    if exc_positions.size:
        exc_arr = np.zeros(exc_positions.size, dtype=_EXC_DTYPE)
        exc_arr["pos"] = exc_positions.astype("<u8")
        exc_arr["ch"] = arr[exc_positions]
        exceptions = exc_arr.tobytes()
        codes[exc_mask] = 0
    else:
        exceptions = b""

    pad = (-n) % 4
    if pad:
        codes = np.concatenate([codes, np.zeros(pad, dtype=np.uint8)])
    codes = codes.reshape(-1, 4)
    packed = (
        codes[:, 0] | (codes[:, 1] << 2) | (codes[:, 2] << 4) | (codes[:, 3] << 6)
    ).astype(np.uint8)
    return packed.tobytes(), exceptions


def _unpack_seq(packed: bytes, exceptions: bytes, n: int) -> bytes:
    packed_arr = np.frombuffer(packed, dtype=np.uint8)
    codes = np.empty(packed_arr.size * 4, dtype=np.uint8)
    codes[0::4] = packed_arr & 0x3
    codes[1::4] = (packed_arr >> 2) & 0x3
    codes[2::4] = (packed_arr >> 4) & 0x3
    codes[3::4] = (packed_arr >> 6) & 0x3
    codes = codes[:n]

    bases = _CODE_LUT[codes].copy()
    if exceptions:
        exc_arr = np.frombuffer(exceptions, dtype=_EXC_DTYPE)
        bases[exc_arr["pos"]] = exc_arr["ch"]
    return bases.tobytes()


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


def compress_file(in_path, out_path, threads: int = 0, try_xz: bool = True, lossy_quality: bool = False):
    ids, plus, seqs, quals, lengths = [], [], [], [], []
    for id_line, seq_line, plus_line, qual_line in _parse_fastq(in_path):
        if len(seq_line) != len(qual_line):
            raise ValueError(f"sequence/quality length mismatch for read {id_line!r}")
        ids.append(id_line)
        plus.append(plus_line)
        seqs.append(seq_line)
        quals.append(qual_line)
        lengths.append(len(seq_line))

    ids_blob = b"\n".join(ids)
    plus_blob = b"\n".join(plus)
    seq_concat = b"".join(seqs)
    qual_concat = b"".join(quals)
    lengths_blob = struct.pack(f"<{len(lengths)}I", *lengths)

    if lossy_quality:
        qual_concat = _bin_quality(qual_concat)

    seq_packed, seq_exceptions = _pack_seq(seq_concat)

    streams_raw = [ids_blob, plus_blob, seq_packed, seq_exceptions, lengths_blob, qual_concat]
    streams = [backend.compress_bytes(s, threads=threads, try_xz=try_xz) for s in streams_raw]

    meta = struct.pack("<BQQ", 1 if lossy_quality else 0, len(seqs), len(seq_concat))
    meta_backend_id, meta_blob = backend.compress_bytes(meta, threads=threads, try_xz=False)

    container.write_container(
        out_path, container.FORMAT_FASTQ, [(meta_backend_id, meta_blob)] + streams
    )

    orig_size = sum(len(x) + 1 for x in ids) * 2 + sum(len(x) for x in seqs) + sum(len(x) for x in quals) + len(seqs) * 2
    return orig_size, sum(len(b) for _, b in [(meta_backend_id, meta_blob)] + streams) + 18


def decompress_file(in_path, out_path):
    format_id, streams = container.read_container(in_path)
    if format_id != container.FORMAT_FASTQ:
        raise ValueError("not a fastq-format .bioz file")

    (meta_backend, meta_blob), *rest = streams
    meta = backend.decompress_bytes(meta_backend, meta_blob)
    lossy_quality, n_reads, total_bases = struct.unpack("<BQQ", meta)

    ids_c, plus_c, seq_packed_c, seq_exc_c, lengths_c, qual_c = rest
    ids_blob = backend.decompress_bytes(*ids_c)
    plus_blob = backend.decompress_bytes(*plus_c)
    seq_packed = backend.decompress_bytes(*seq_packed_c)
    seq_exceptions = backend.decompress_bytes(*seq_exc_c)
    lengths_blob = backend.decompress_bytes(*lengths_c)
    qual_concat = backend.decompress_bytes(*qual_c)

    lengths = struct.unpack(f"<{n_reads}I", lengths_blob)
    seq_concat = _unpack_seq(seq_packed, seq_exceptions, total_bases)

    ids = ids_blob.split(b"\n") if n_reads else []
    plus = plus_blob.split(b"\n") if n_reads else []

    with open(out_path, "wb") as f:
        offset = 0
        for i in range(n_reads):
            n = lengths[i]
            seq = seq_concat[offset : offset + n]
            qual = qual_concat[offset : offset + n]
            offset += n
            f.write(b"@" + ids[i] + b"\n")
            f.write(seq + b"\n")
            f.write(b"+" + plus[i] + b"\n")
            f.write(qual + b"\n")

    return lossy_quality
