"""Fallback codec for any file that isn't a recognized bioinformatics format:
whole file goes through the best-of backend chooser as a single stream."""

from . import backend
from . import container


def compress_file(in_path, out_path, threads: int = 0, try_xz: bool = True):
    with open(in_path, "rb") as f:
        data = f.read()
    backend_id, blob = backend.compress_bytes(data, threads=threads, try_xz=try_xz)
    container.write_container(out_path, container.FORMAT_GENERIC, [(backend_id, blob)])
    return len(data), len(blob) + 18  # +approx container overhead


def decompress_file(in_path, out_path):
    format_id, streams = container.read_container(in_path)
    if format_id != container.FORMAT_GENERIC:
        raise ValueError("not a generic-format .bioz file")
    (backend_id, blob), = streams
    data = backend.decompress_bytes(backend_id, blob)
    with open(out_path, "wb") as f:
        f.write(data)
