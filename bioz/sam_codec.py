"""
Columnar SAM/BAM codec -- no reference genome required.

Splits SAM alignment records into per-field streams (QNAME, FLAG, RNAME,
POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, optional-tags) instead of
compressing the interleaved tab-separated text/binary directly. Each column
is far more self-similar on its own (e.g. all the CIGAR strings together,
all the RNAME values together) than the mixed per-record stream is.

This is the fallback for aligned data when no reference is available for
CRAM (bam_codec.py); CRAM (reference-based) beats this when you have the
reference, since it can avoid storing redundant sequence entirely.

Correctness note: .sam input round-trips byte-exact. .bam input is
regenerated via `samtools view -b` from the reconstructed SAM text --
record-for-record identical (verified in tests), but not guaranteed
byte-identical BAM *binary*, since samtools re-encodes it rather than the
original bytes being reproduced bit-for-bit.
"""

import subprocess
from pathlib import Path

from . import backend
from . import container

N_FIELDS = 11  # QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL


def _to_sam_text(in_path) -> bytes:
    suffix = str(in_path).lower()
    if suffix.endswith(".bam") or suffix.endswith(".cram"):
        proc = subprocess.run(["samtools", "view", "-h", str(in_path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if proc.returncode != 0:
            raise RuntimeError(f"samtools view failed: {proc.stderr.decode(errors='replace')}")
        return proc.stdout
    return Path(in_path).read_bytes()


def compress_file(in_path, out_path, threads: int = 0, try_xz: bool = True):
    sam_bytes = _to_sam_text(in_path)
    lines = sam_bytes.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()

    header_lines = []
    i = 0
    while i < len(lines) and lines[i].startswith(b"@"):
        header_lines.append(lines[i])
        i += 1

    fields = [[] for _ in range(N_FIELDS)]
    tags = []
    for line in lines[i:]:
        parts = line.split(b"\t", N_FIELDS)
        if len(parts) < N_FIELDS:
            raise ValueError(f"malformed SAM record (fewer than {N_FIELDS} fields): {line!r}")
        for k in range(N_FIELDS):
            fields[k].append(parts[k])
        tags.append(parts[N_FIELDS] if len(parts) > N_FIELDS else b"")

    header_blob = b"\n".join(header_lines)
    field_blobs = [b"\n".join(f) for f in fields]
    tags_blob = b"\n".join(tags)

    all_blobs = [header_blob] + field_blobs + [tags_blob]
    streams = [backend.compress_bytes(b, threads=threads, try_xz=try_xz) for b in all_blobs]

    n_records = len(lines) - i
    meta = n_records.to_bytes(8, "little") + len(header_lines).to_bytes(8, "little")
    meta_id, meta_blob = backend.compress_bytes(meta, threads=threads, try_xz=False)

    container.write_container(out_path, container.FORMAT_SAM, [(meta_id, meta_blob)] + streams)
    orig = Path(in_path).stat().st_size
    return orig, Path(out_path).stat().st_size


def decompress_file(in_path, out_path, threads: int = 0):
    format_id, streams = container.read_container(in_path)
    if format_id != container.FORMAT_SAM:
        raise ValueError("not a SAM-format .bioz file")

    (meta_id, meta_blob), *rest = streams
    meta = backend.decompress_bytes(meta_id, meta_blob)
    n_records = int.from_bytes(meta[0:8], "little")
    n_header = int.from_bytes(meta[8:16], "little")

    header_c, *field_and_tag_c = rest
    field_c = field_and_tag_c[:N_FIELDS]
    (tags_c,) = field_and_tag_c[N_FIELDS:]

    header_blob = backend.decompress_bytes(*header_c)
    field_blobs = [backend.decompress_bytes(*c) for c in field_c]
    tags_blob = backend.decompress_bytes(*tags_c)

    header_lines = header_blob.split(b"\n") if n_header else []
    fields = [fb.split(b"\n") if n_records else [] for fb in field_blobs]
    tags = tags_blob.split(b"\n") if n_records else []

    sam_lines = list(header_lines)
    for i in range(n_records):
        rec_fields = [fields[k][i] for k in range(N_FIELDS)]
        tag = tags[i]
        if tag:
            rec_fields.append(tag)
        sam_lines.append(b"\t".join(rec_fields))
    sam_bytes = b"\n".join(sam_lines) + (b"\n" if sam_lines else b"")

    out_ext = Path(out_path).suffix.lower()
    if out_ext == ".sam":
        Path(out_path).write_bytes(sam_bytes)
    else:
        proc = subprocess.run(
            ["samtools", "view", "-b", "-o", str(out_path), "-"],
            input=sam_bytes, stderr=subprocess.PIPE,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"samtools SAM->BAM failed: {proc.stderr.decode(errors='replace')}")
