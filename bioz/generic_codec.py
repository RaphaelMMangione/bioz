"""Fallback codec for any file that isn't a recognized bioinformatics format:
whole file goes through the best-of backend chooser as a single stream,
file-to-file so peak memory doesn't scale with the file's size."""

from pathlib import Path

from . import backend
from . import container


def compress_file(in_path, out_path, threads: int = 0, try_xz: bool = True, ultra: bool = False):
    in_path = Path(in_path)
    backend_id, tmp_path, comp_size = backend.compress_file(in_path, threads=threads, try_xz=try_xz, ultra=ultra)
    try:
        container.write_container(out_path, container.FORMAT_GENERIC, [(backend_id, tmp_path)])
    finally:
        tmp_path.unlink(missing_ok=True)
    return in_path.stat().st_size, comp_size + 18  # +approx container overhead


def decompress_file(in_path, out_path):
    format_id, index = container.iter_container_streams(in_path)
    if format_id != container.FORMAT_GENERIC:
        raise ValueError("not a generic-format .bioz file")
    (backend_id, offset, length), = index

    import tempfile
    tmp_path = Path(tempfile.mkstemp(prefix="bioz_gen_extract_")[1])
    try:
        container.extract_stream_to_file(in_path, offset, length, tmp_path)
        backend.decompress_to_file(backend_id, tmp_path, out_path)
    finally:
        tmp_path.unlink(missing_ok=True)
