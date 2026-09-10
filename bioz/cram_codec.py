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

samtools always reads/writes real files here (never through a Python
`bytes` object) so peak memory doesn't scale with alignment file size.
"""

import hashlib
import subprocess
import tempfile
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


def compress_file(in_path, out_path, ref_path, threads: int = 0, try_xz: bool = True, ultra: bool = False):
    ref_path = Path(ref_path)
    if not ref_path.exists():
        raise FileNotFoundError(f"reference not found: {ref_path}")

    nthreads = max(1, threads or 1)
    cram_tmp = Path(tempfile.mkstemp(prefix="bioz_cram_")[1])
    try:
        proc = subprocess.run(
            ["samtools", "view", "-C", "-T", str(ref_path), "-@", str(nthreads),
             "-o", str(cram_tmp), str(in_path)],
            stderr=subprocess.PIPE,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"samtools BAM/SAM->CRAM conversion failed: {proc.stderr.decode(errors='replace')}")

        data_backend, data_tmp, data_size = backend.compress_file(cram_tmp, threads=threads, try_xz=try_xz, ultra=ultra)
        try:
            ref_abspath = str(ref_path.resolve()).encode()
            ref_md5 = _file_md5(ref_path).encode()
            meta = len(ref_abspath).to_bytes(4, "little") + ref_abspath + ref_md5
            meta_backend, meta_blob = backend.compress_bytes(meta, threads=threads, try_xz=False)

            container.write_container(
                out_path, container.FORMAT_CRAM, [(meta_backend, meta_blob), (data_backend, data_tmp)]
            )
        finally:
            data_tmp.unlink(missing_ok=True)
    finally:
        cram_tmp.unlink(missing_ok=True)

    orig = Path(in_path).stat().st_size
    return orig, Path(out_path).stat().st_size


def decompress_file(in_path, out_path, ref_path=None, threads: int = 0):
    format_id, index = container.iter_container_streams(in_path)
    if format_id != container.FORMAT_CRAM:
        raise ValueError("not a CRAM-format .bioz file")

    (meta_backend, meta_off, meta_len), (data_backend, data_off, data_len) = index
    meta_tmp = Path(tempfile.mkstemp(prefix="bioz_cram_meta_")[1])
    try:
        container.extract_stream_to_file(in_path, meta_off, meta_len, meta_tmp)
        meta = meta_tmp.read_bytes()
    finally:
        meta_tmp.unlink(missing_ok=True)
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

    data_extract_tmp = Path(tempfile.mkstemp(prefix="bioz_cram_data_")[1])
    cram_tmp = Path(tempfile.mkstemp(prefix="bioz_cram_out_")[1])
    try:
        container.extract_stream_to_file(in_path, data_off, data_len, data_extract_tmp)
        backend.decompress_to_file(data_backend, data_extract_tmp, cram_tmp)

        out_ext = Path(out_path).suffix.lower()
        view_flag = "-h" if out_ext == ".sam" else "-b"
        nthreads = max(1, threads or 1)
        proc = subprocess.run(
            ["samtools", "view", view_flag, "-T", str(ref), "-@", str(nthreads),
             "-o", str(out_path), str(cram_tmp)],
            stderr=subprocess.PIPE,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"samtools CRAM->BAM/SAM conversion failed: {proc.stderr.decode(errors='replace')}")
    finally:
        data_extract_tmp.unlink(missing_ok=True)
        cram_tmp.unlink(missing_ok=True)
