"""
Generic byte-stream compression backend. Not a novel entropy coder -- this
picks the best of the strongest available general-purpose compressors
(zstd --ultra -22 --long, xz -9e) per blob, so higher layers (fastq codec,
generic file codec) don't have to know or care which one won.
"""

import shutil
import subprocess

ZSTD = shutil.which("zstd")
XZ = shutil.which("xz")

# backend id bytes stored in the container so decompress knows what to invoke
BACKEND_ZSTD = 1
BACKEND_XZ = 2
BACKEND_STORE = 0  # no compression helped -- stored raw


def _run(cmd, data: bytes) -> bytes:
    proc = subprocess.run(cmd, input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {proc.stderr.decode(errors='replace')}")
    return proc.stdout


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
    to storing raw if neither helps (e.g. already-compressed input)."""
    if threads <= 0:
        import os
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
