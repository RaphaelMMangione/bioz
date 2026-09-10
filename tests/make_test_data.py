"""Generates small synthetic test files under tests/test_data/ -- nothing here
touches any real user data."""
import gzip
import random
from pathlib import Path

random.seed(42)
OUT = Path(__file__).parent / "test_data"
OUT.mkdir(exist_ok=True)

BASES = "ACGT"


def make_fastq(path, n_reads, len_range, n_rate=0.0, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "wt") as f:
        for i in range(n_reads):
            L = random.randint(*len_range)
            seq = []
            for _ in range(L):
                if random.random() < n_rate:
                    seq.append("N")
                else:
                    seq.append(random.choice(BASES))
            seq = "".join(seq)
            qual = "".join(chr(33 + random.randint(2, 40)) for _ in range(L))
            f.write(f"@read{i}-{random.randint(0,10**8)} runid=abc123 ch={i%512} start_time=2026-09-01T00:00:00Z\n")
            f.write(seq + "\n")
            f.write("+\n")
            f.write(qual + "\n")


# small: quick correctness/roundtrip checks
make_fastq(OUT / "small.fastq", n_reads=200, len_range=(50, 300), n_rate=0.01)

# medium: short-read-like (mirrors the real dataset's ~150-300bp reads), for a realistic benchmark
make_fastq(OUT / "medium_shortreads.fastq", n_reads=20000, len_range=(80, 400), n_rate=0.005)

# a plain generic (non-bioinformatics) file: repetitive + some random text, to exercise the generic path
with open(OUT / "generic.txt", "w") as f:
    for _ in range(20000):
        f.write("the quick brown fox jumps over the lazy dog " * random.randint(1, 3))
        f.write(str(random.random()) + "\n")

print("test data written to", OUT)
for p in sorted(OUT.iterdir()):
    print(" ", p.name, p.stat().st_size, "bytes")
