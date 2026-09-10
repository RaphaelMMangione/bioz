"""
Generic byte-stream compression backend. Not a novel entropy coder -- this
picks the best of the strongest available general-purpose compressors
(zstd --ultra -22 --long, xz -9e) per stream, so higher layers (fastq
codec, generic file codec) don't have to know or care which one won.

Everything here is file-to-file (never holds a whole stream as a single
Python `bytes` object) so peak memory stays bounded regardless of input
size -- the OS pipes bytes from the source file straight into the
compressor subprocess and from the subprocess straight into the output
file. `compress_bytes`/`decompress_bytes` remain for small, bounded blobs
(metadata headers, etc.) where an in-memory bytes object is fine.
"""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

ZSTD = shutil.which("zstd")
XZ = shutil.which("xz")

# backend id bytes stored in the container so decompress knows what to invoke
BACKEND_ZSTD = 1
BACKEND_XZ = 2
BACKEND_STORE = 0  # no compression helped -- stored raw

_CHUNK = 1 << 20  # 1MiB, used for any manual streaming copy in this module


def _run(cmd, data: bytes) -> bytes:
    proc = subprocess.run(cmd, input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {proc.stderr.decode(errors='replace')}")
    return proc.stdout


def _run_file_to_file(cmd, in_path, out_path):
    with open(in_path, "rb") as fin, open(out_path, "wb") as fout:
        proc = subprocess.run(cmd, stdin=fin, stdout=fout, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {proc.stderr.decode(errors='replace')}")


def _zstd_compress(data: bytes, threads: int) -> bytes:
    return _run([ZSTD, "-q", "--ultra", "-22", "--long=27", f"-T{threads}", "-c"], data)


def _zstd_decompress(data: bytes) -> bytes:
    return _run([ZSTD, "-q", "-d", "--long=27", "-c"], data)


def _xz_compress(data: bytes, threads: int) -> bytes:
    return _run([XZ, "-q", "-9e", f"-T{threads}", "-c"], data)


def _xz_decompress(data: bytes) -> bytes:
    return _run([XZ, "-q", "-d", "-c"], data)


def compress_bytes(data: bytes, threads: int = 0, try_xz: bool = True) -> tuple:
    """Compress `data`, returning (backend_id: int, compressed: bytes).
    Picks whichever backend actually produces the smaller output; falls back
    to storing raw if neither helps (e.g. already-compressed input). For
    small, bounded blobs only (metadata) -- see compress_file for anything
    that could be large."""
    if threads <= 0:
        threads = max(1, os.cpu_count() or 1)

    candidates = [(BACKEND_STORE, data)]

    if ZSTD and len(data) > 0:
        try:
            candidates.append((BACKEND_ZSTD, _zstd_compress(data, threads)))
        except Exception:
            pass

    if try_xz and XZ and len(data) > 0:
        try:
            candidates.append((BACKEND_XZ, _xz_compress(data, threads)))
        except Exception:
            pass

    backend, out = min(candidates, key=lambda c: len(c[1]))
    return backend, out


def decompress_bytes(backend: int, data: bytes) -> bytes:
    if backend == BACKEND_STORE:
        return data
    elif backend == BACKEND_ZSTD:
        return _zstd_decompress(data)
    elif backend == BACKEND_XZ:
        return _xz_decompress(data)
    else:
        raise ValueError(f"unknown backend id {backend}")


def compress_file(in_path, threads: int = 0, try_xz: bool = True, tmp_dir=None) -> tuple:
    """Compress the file at `in_path`, returning (backend_id, out_path,
    size). Runs each candidate backend reading/writing files directly (no
    Python-side buffering of the stream), so peak memory is independent of
    `in_path`'s size. Caller is responsible for deleting `out_path` once
    its contents have been consumed (e.g. copied into a container)."""
    if threads <= 0:
        threads = max(1, os.cpu_count() or 1)
    in_path = Path(in_path)
    in_size = in_path.stat().st_size

    candidates = [(BACKEND_STORE, in_path, in_size)]  # store: reuse input, don't copy yet

    if ZSTD and in_size > 0:
        out = Path(tempfile.mkstemp(prefix="bioz_zstd_", dir=tmp_dir)[1])
        try:
            _run_file_to_file([ZSTD, "-q", "--ultra", "-22", "--long=27", f"-T{threads}", "-c"], in_path, out)
            candidates.append((BACKEND_ZSTD, out, out.stat().st_size))
        except Exception:
            out.unlink(missing_ok=True)

    if try_xz and XZ and in_size > 0:
        out = Path(tempfile.mkstemp(prefix="bioz_xz_", dir=tmp_dir)[1])
        try:
            _run_file_to_file([XZ, "-q", "-9e", f"-T{threads}", "-c"], in_path, out)
            candidates.append((BACKEND_XZ, out, out.stat().st_size))
        except Exception:
            out.unlink(missing_ok=True)

    winner_backend, winner_path, winner_size = min(candidates, key=lambda c: c[2])

    # clean up the losers (but never delete the original input file)
    for cand_backend, cand_path, _ in candidates:
        if cand_path != winner_path and cand_path != in_path:
            Path(cand_path).unlink(missing_ok=True)

    if winner_path == in_path:
        # STORE won: copy to a real temp file so the caller can always
        # treat the result as "ours to delete" uniformly.
        store_path = Path(tempfile.mkstemp(prefix="bioz_store_", dir=tmp_dir)[1])
        with open(in_path, "rb") as fin, open(store_path, "wb") as fout:
            shutil.copyfileobj(fin, fout, _CHUNK)
        winner_path = store_path

    return winner_backend, winner_path, winner_size


def decompress_to_file(backend: int, in_path, out_path):
    """Decompress the file at `in_path` (compressed with the given backend
    id) directly to `out_path`, file-to-file."""
    if backend == BACKEND_STORE:
        with open(in_path, "rb") as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout, _CHUNK)
    elif backend == BACKEND_ZSTD:
        _run_file_to_file([ZSTD, "-q", "-d", "--long=27", "-c"], in_path, out_path)
    elif backend == BACKEND_XZ:
        _run_file_to_file([XZ, "-q", "-d", "-c"], in_path, out_path)
    else:
        raise ValueError(f"unknown backend id {backend}")
