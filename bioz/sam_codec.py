"""
Columnar SAM/BAM codec -- no reference genome required.

Splits SAM alignment records into per-field streams (QNAME, FLAG, RNAME,
POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, optional-tags) instead of
compressing the interleaved tab-separated text/binary directly. Each column
is far more self-similar on its own (e.g. all the CIGAR strings together,
all the RNAME values together) than the mixed per-record stream is.

This is the fallback for aligned data when no reference is available for
CRAM (cram_codec.py); CRAM (reference-based) beats this when you have the
reference, since it can avoid storing redundant sequence entirely.

Correctness note: .sam input round-trips byte-exact. .bam input is
regenerated via `samtools view -b` from the reconstructed SAM text --
record-for-record identical (verified in tests), but not guaranteed
byte-identical BAM *binary*, since samtools re-encodes it rather than the
original bytes being reproduced bit-for-bit.

Both directions process records in bounded-size chunks and stream every
per-field stream through temp files, so peak memory stays roughly constant
regardless of the alignment file's size (samtools itself reads/writes real
files too, never a Python-buffered blob).
"""

import struct
import subprocess
import tempfile
from pathlib import Path

from . import backend
from . import container

N_FIELDS = 11  # QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
CHUNK_RECORDS = 200_000


def _sam_lines(in_path):
    """Yields raw SAM lines (bytes, no trailing newline). For .bam/.cram,
    samtools writes to a temp file first (not captured into a Python bytes
    object), which is then read line-by-line and cleaned up."""
    suffix = str(in_path).lower()
    if suffix.endswith(".bam") or suffix.endswith(".cram"):
        tmp = Path(tempfile.mkstemp(prefix="bioz_sam_view_")[1])
        try:
            # redirect at the OS level (not samtools' own -o) so its @PG CL:
            # tag records the same invocation as capturing stdout would --
            # streams straight to disk either way, never through Python.
            with open(tmp, "wb") as tmp_f:
                proc = subprocess.run(["samtools", "view", "-h", str(in_path)],
                                       stdout=tmp_f, stderr=subprocess.PIPE)
            if proc.returncode != 0:
                raise RuntimeError(f"samtools view failed: {proc.stderr.decode(errors='replace')}")
            with open(tmp, "rb") as f:
                for line in f:
                    yield line.rstrip(b"\n")
        finally:
            tmp.unlink(missing_ok=True)
    else:
        with open(in_path, "rb") as f:
            for line in f:
                yield line.rstrip(b"\n")


def compress_file(in_path, out_path, threads: int = 0, try_xz: bool = True, ultra: bool = False):
    tmp_dir = Path(tempfile.mkdtemp(prefix="bioz_sam_"))
    field_names = [f"field{k}" for k in range(N_FIELDS)]
    paths = {name: tmp_dir / name for name in field_names + ["header", "tags"]}
    n_records = 0
    n_header = 0

    try:
        files = {name: open(p, "wb") for name, p in paths.items()}
        first_field = {name: True for name in field_names + ["tags"]}
        first_header = True
        try:
            for line in _sam_lines(in_path):
                if line.startswith(b"@"):
                    files["header"].write((b"" if first_header else b"\n") + line)
                    first_header = False
                    n_header += 1
                    continue
                parts = line.split(b"\t", N_FIELDS)
                if len(parts) < N_FIELDS:
                    raise ValueError(f"malformed SAM record (fewer than {N_FIELDS} fields): {line!r}")
                for k, name in enumerate(field_names):
                    files[name].write((b"" if first_field[name] else b"\n") + parts[k])
                    first_field[name] = False
                tag = parts[N_FIELDS] if len(parts) > N_FIELDS else b""
                files["tags"].write((b"" if first_field["tags"] else b"\n") + tag)
                first_field["tags"] = False
                n_records += 1
        finally:
            for f in files.values():
                f.close()

        stream_order = ["header"] + field_names + ["tags"]
        streams = []
        for name in stream_order:
            backend_id, comp_path, _size = backend.compress_file(paths[name], threads=threads, try_xz=try_xz, ultra=ultra)
            streams.append((backend_id, comp_path))

        meta = n_records.to_bytes(8, "little") + n_header.to_bytes(8, "little")
        meta_id, meta_blob = backend.compress_bytes(meta, threads=threads, try_xz=False)

        container.write_container(out_path, container.FORMAT_SAM, [(meta_id, meta_blob)] + streams)
        for _, comp_path in streams:
            Path(comp_path).unlink(missing_ok=True)

        orig = Path(in_path).stat().st_size
        return orig, Path(out_path).stat().st_size
    finally:
        for p in paths.values():
            p.unlink(missing_ok=True)
        tmp_dir.rmdir()


def decompress_file(in_path, out_path, threads: int = 0):
    format_id, index = container.iter_container_streams(in_path)
    if format_id != container.FORMAT_SAM:
        raise ValueError("not a SAM-format .bioz file")

    (meta_id, meta_off, meta_len), *rest_idx = index
    tmp_dir = Path(tempfile.mkdtemp(prefix="bioz_sam_dec_"))
    try:
        meta_c = tmp_dir / "meta.c"
        container.extract_stream_to_file(in_path, meta_off, meta_len, meta_c)
        meta = backend.decompress_bytes(meta_id, meta_c.read_bytes())
        meta_c.unlink(missing_ok=True)
        n_records = int.from_bytes(meta[0:8], "little")
        n_header = int.from_bytes(meta[8:16], "little")

        field_names = [f"field{k}" for k in range(N_FIELDS)]
        stream_order = ["header"] + field_names + ["tags"]
        decompressed = {}
        for name, (b_id, off, length) in zip(stream_order, rest_idx):
            ex_path = tmp_dir / f"{name}.c"
            container.extract_stream_to_file(in_path, off, length, ex_path)
            dec_path = tmp_dir / name
            backend.decompress_to_file(b_id, ex_path, dec_path)
            ex_path.unlink(missing_ok=True)
            decompressed[name] = dec_path

        out_ext = Path(out_path).suffix.lower()
        sam_target = Path(out_path) if out_ext == ".sam" else tmp_dir / "reconstructed.sam"

        field_files = [open(decompressed[name], "rb") for name in field_names]
        tags_f = open(decompressed["tags"], "rb")
        try:
            with open(sam_target, "wb") as out_f, open(decompressed["header"], "rb") as hdr_f:
                for i, line in enumerate(hdr_f):
                    out_f.write(line if line.endswith(b"\n") else line + b"\n")
                for _ in range(n_records):
                    rec_fields = [f.readline().rstrip(b"\n") for f in field_files]
                    tag = tags_f.readline().rstrip(b"\n")
                    if tag:
                        rec_fields.append(tag)
                    out_f.write(b"\t".join(rec_fields) + b"\n")
        finally:
            for f in field_files:
                f.close()
            tags_f.close()

        if out_ext != ".sam":
            proc = subprocess.run(
                ["samtools", "view", "-b", "-o", str(out_path), str(sam_target)],
                stderr=subprocess.PIPE,
            )
            if proc.returncode != 0:
                raise RuntimeError(f"samtools SAM->BAM failed: {proc.stderr.decode(errors='replace')}")
    finally:
        for p in tmp_dir.iterdir():
            p.unlink(missing_ok=True)
        tmp_dir.rmdir()
