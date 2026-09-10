"""
Minimal container format for .bioz archives.

Layout:
  magic       b"BIOZ1"          (5 bytes)
  format_id   1 byte            (0 = generic passthrough, 1 = fastq-aware)
  n_streams   uint32 LE
  for each stream:
      backend_id  1 byte
      length      uint64 LE
      bytes       <length>
"""

import shutil
import struct
from pathlib import Path

MAGIC = b"BIOZ1"

FORMAT_GENERIC = 0
FORMAT_FASTQ = 1
FORMAT_SAM = 2
FORMAT_CRAM = 3
FORMAT_SIGNAL = 4

_CHUNK = 1 << 20  # 1MiB


def write_container(path, format_id: int, streams: list) -> None:
    """streams: list of (backend_id, blob) tuples, in a fixed order agreed
    with the reader. `blob` may be `bytes` (kept for small metadata streams)
    or a `Path`/path-string, whose file contents are streamed straight into
    the container in chunks -- the whole stream is never held in memory at
    once, so container size can exceed available RAM."""
    with open(path, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<B", format_id))
        f.write(struct.pack("<I", len(streams)))
        for backend_id, blob in streams:
            if isinstance(blob, (bytes, bytearray)):
                f.write(struct.pack("<B", backend_id))
                f.write(struct.pack("<Q", len(blob)))
                f.write(blob)
            else:
                blob_path = Path(blob)
                length = blob_path.stat().st_size
                f.write(struct.pack("<B", backend_id))
                f.write(struct.pack("<Q", length))
                with open(blob_path, "rb") as src:
                    shutil.copyfileobj(src, f, _CHUNK)


def read_container(path):
    """Reads every stream fully into memory -- fine for small containers
    (metadata-only formats) or when the caller already needs everything at
    once. For large formats, use iter_container_streams + extract_stream_to_file
    instead so no single stream needs to fit in RAM."""
    with open(path, "rb") as f:
        magic = f.read(5)
        if magic != MAGIC:
            raise ValueError(f"{path} is not a .bioz file (bad magic)")
        (format_id,) = struct.unpack("<B", f.read(1))
        (n_streams,) = struct.unpack("<I", f.read(4))
        streams = []
        for _ in range(n_streams):
            (backend_id,) = struct.unpack("<B", f.read(1))
            (length,) = struct.unpack("<Q", f.read(8))
            blob = f.read(length)
            streams.append((backend_id, blob))
        return format_id, streams


def iter_container_streams(path):
    """Reads only the header/index, not the stream contents -- returns
    (format_id, [(backend_id, offset, length), ...]). Use extract_stream_to_file
    to pull out any one stream's raw (still-compressed) bytes without
    loading the others, or loading this one fully into memory either."""
    with open(path, "rb") as f:
        magic = f.read(5)
        if magic != MAGIC:
            raise ValueError(f"{path} is not a .bioz file (bad magic)")
        (format_id,) = struct.unpack("<B", f.read(1))
        (n_streams,) = struct.unpack("<I", f.read(4))
        index = []
        offset = f.tell()
        for _ in range(n_streams):
            f.seek(offset)
            backend_id_b = f.read(1)
            length_b = f.read(8)
            (backend_id,) = struct.unpack("<B", backend_id_b)
            (length,) = struct.unpack("<Q", length_b)
            data_offset = offset + 1 + 8
            index.append((backend_id, data_offset, length))
            offset = data_offset + length
        return format_id, index


def extract_stream_to_file(container_path, offset, length, out_path):
    """Chunked copy of container_path[offset:offset+length] to out_path,
    without ever holding more than one chunk in memory."""
    with open(container_path, "rb") as src, open(out_path, "wb") as dst:
        src.seek(offset)
        remaining = length
        while remaining > 0:
            chunk = src.read(min(_CHUNK, remaining))
            if not chunk:
                raise ValueError(f"{container_path}: truncated stream (expected {length} bytes at offset {offset})")
            dst.write(chunk)
            remaining -= len(chunk)
