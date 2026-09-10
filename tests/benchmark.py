"""Benchmarks bioz against gzip -9, zstd -19, xz -9 on the test data.
Prints a plain comparison table -- no claims, just measured numbers."""
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).parent
DATA = HERE / "test_data"
WORK = HERE / "bench_work"
WORK.mkdir(exist_ok=True)

FILES = ["small.fastq", "medium_shortreads.fastq", "real_slice.fastq", "generic.txt"]


def timed(cmd):
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True)
    dt = time.time() - t0
    if r.returncode != 0:
        raise RuntimeError(f"{cmd} failed: {r.stderr.decode(errors='replace')}")
    return dt


def size(p):
    return Path(p).stat().st_size


rows = []
for fname in FILES:
    src = DATA / fname
    orig = size(src)
    row = {"file": fname, "orig": orig}

    gz_out = WORK / (fname + ".gz")
    dt = timed(["bash", "-c", f"gzip -9 -c '{src}' > '{gz_out}'"])
    row["gzip9"] = (size(gz_out), dt)

    zstd_out = WORK / (fname + ".zst")
    dt = timed(["bash", "-c", f"zstd -q -19 -T0 -c '{src}' > '{zstd_out}'"])
    row["zstd19"] = (size(zstd_out), dt)

    xz_out = WORK / (fname + ".xz")
    dt = timed(["bash", "-c", f"xz -q -9e -T0 -c '{src}' > '{xz_out}'"])
    row["xz9e"] = (size(xz_out), dt)

    bioz_out = WORK / (fname + ".bioz")
    dt = timed(["bioz", "compress", str(src), "-o", str(bioz_out)])
    row["bioz"] = (size(bioz_out), dt)

    bioz_lossy_out = WORK / (fname + ".lossy.bioz")
    dt = timed(["bioz", "compress", str(src), "-o", str(bioz_lossy_out), "--lossy-quality"])
    row["bioz_lossy"] = (size(bioz_lossy_out), dt)

    rows.append(row)

print(f"{'file':28} {'orig':>10} {'gzip -9':>16} {'zstd -19':>16} {'xz -9e':>16} {'bioz':>16} {'bioz --lossy-quality':>22}")
for row in rows:
    def fmt(key):
        s, dt = row[key]
        ratio = row["orig"] / s if s else float("inf")
        return f"{s:>9} ({ratio:4.2f}x,{dt:5.1f}s)"

    print(f"{row['file']:28} {row['orig']:>10} {fmt('gzip9'):>16} {fmt('zstd19'):>16} {fmt('xz9e'):>16} {fmt('bioz'):>16} {fmt('bioz_lossy'):>22}")
