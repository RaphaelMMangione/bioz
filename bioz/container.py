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

import struct

MAGIC = b"BIOZ1"

FORMAT_GENERIC = 0
FORMAT_FASTQ = 1
FORMAT_SAM = 2
FORMAT_CRAM = 3
FORMAT_SIGNAL = 4


def write_container(path, format_id: int, streams: list) -> None:
    """streams: list of (backend_id, compressed_bytes) tuples, in a fixed
    order agreed with the reader (e.g. for fastq: ids, plus, seq_packed,
    seq_exceptions, seq_lengths, qual)."""
    with open(path, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<B", format_id))
        f.write(struct.pack("<I", len(streams)))
        for backend_id, blob in streams:
            f.write(struct.pack("<B", backend_id))
            f.write(struct.pack("<Q", len(blob)))
            f.write(blob)


def read_container(path):
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
