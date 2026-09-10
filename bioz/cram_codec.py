"""
Reference-based codec for aligned reads (SAM/BAM -> CRAM) -- requires a
reference FASTA. This is the actual proven best-in-class technique for
aligned data: since every read's sequence is already implied by where it
maps on a known reference plus its differences (SNPs/indels), CRAM can
avoid storing most of the sequence at all. Nothing generic (gzip/zstd/xz)
can do this, because they don't know the reference exists.

The reference's absolute path and MD5 are stored in the .bioz container so
decompression can verify the same reference is being used -- CRAM decoding
with the wrong reference would silently reconstruct the wrong sequence
rather than erroring, so this is checked rather than trusted.
"""

import hashlib
import subprocess
from pathlib import Path

from . import backend
from . import container


def _file_md5(path, chunk_size=1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def compress_file(in_path, out_path, ref_path, threads: int = 0, try_xz: bool = True):
    ref_path = Path(ref_path)
    if not ref_path.exists():
        raise FileNotFoundError(f"reference not found: {ref_path}")

    nthreads = max(1, threads or 1)
    proc = subprocess.run(
        ["samtools", "view", "-C", "-T", str(ref_path), "-@", str(nthreads), str(in_path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"samtools BAM/SAM->CRAM conversion failed: {proc.stderr.decode(errors='replace')}")
    cram_bytes = proc.stdout

    data_backend, data_blob = backend.compress_bytes(cram_bytes, threads=threads, try_xz=try_xz)

    ref_abspath = str(ref_path.resolve()).encode()
    ref_md5 = _file_md5(ref_path).encode()
    meta = len(ref_abspath).to_bytes(4, "little") + ref_abspath + ref_md5
    meta_backend, meta_blob = backend.compress_bytes(meta, threads=threads, try_xz=False)

    container.write_container(
        out_path, container.FORMAT_CRAM, [(meta_backend, meta_blob), (data_backend, data_blob)]
    )
    orig = Path(in_path).stat().st_size
    return orig, Path(out_path).stat().st_size


def decompress_file(in_path, out_path, ref_path=None, threads: int = 0):
    format_id, streams = container.read_container(in_path)
    if format_id != container.FORMAT_CRAM:
        raise ValueError("not a CRAM-format .bioz file")

    (meta_backend, meta_blob), (data_backend, data_blob) = streams
    meta = backend.decompress_bytes(meta_backend, meta_blob)
    n = int.from_bytes(meta[:4], "little")
    stored_ref_path = meta[4 : 4 + n].decode()
    stored_ref_md5 = meta[4 + n :].decode()

    ref = Path(ref_path) if ref_path else Path(stored_ref_path)
    if not ref.exists():
        raise FileNotFoundError(
            f"reference used at compress time is not available: {stored_ref_path}\n"
            f"pass --ref /path/to/the/same/reference.fa to decompress"
        )
    if _file_md5(ref) != stored_ref_md5:
        raise ValueError(
            f"reference at {ref} does not match the one used to compress this file (MD5 mismatch) -- "
            f"decoding CRAM against the wrong reference would silently reconstruct the wrong sequence, refusing"
        )

    cram_bytes = backend.decompress_bytes(data_backend, data_blob)

    out_ext = Path(out_path).suffix.lower()
    view_flag = "-h" if out_ext == ".sam" else "-b"
    nthreads = max(1, threads or 1)
    proc = subprocess.run(
        ["samtools", "view", view_flag, "-T", str(ref), "-@", str(nthreads), "-o", str(out_path), "-"],
        input=cram_bytes, stderr=subprocess.PIPE,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"samtools CRAM->BAM/SAM conversion failed: {proc.stderr.decode(errors='replace')}")
