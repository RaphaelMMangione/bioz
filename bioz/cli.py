import argparse
import gzip
import sys
import time
from pathlib import Path

from . import fastq_codec, generic_codec, sam_codec, cram_codec, container


def _detect_format(path: Path, force_generic: bool) -> str:
    """Returns one of: 'fastq', 'sam', 'generic'."""
    if force_generic:
        return "generic"
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith((".fastq", ".fq")):
        return "fastq"
    if name.endswith((".bam", ".sam", ".cram")):
        return "sam"
    if name.endswith((".fasta", ".fa")):
        return "generic"
    # unknown extension: peek at the first bytes to guess
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rb") as f:
            first = f.read(4)
        if first[:1] == b"@" and path.suffix != ".sam":
            return "fastq"
        if first == b"BAM\x01" or first[:4] == b"CRAM":
            return "sam"
    except Exception:
        pass
    return "generic"


def _human(n: int) -> str:
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}PB"


def _compress_one(in_path: Path, out_path: Path, args) -> str:
    """Compresses a single file, dispatching on its detected format.
    Returns the mode string used (for reporting)."""
    fmt = _detect_format(in_path, args.generic)

    if fmt == "fastq":
        fastq_codec.compress_file(
            in_path, out_path, threads=args.threads, try_xz=args.max, lossy_quality=args.lossy_quality
        )
        return "fastq-aware" + (", lossy quality" if args.lossy_quality else ", lossless")
    elif fmt == "sam" and args.ref:
        cram_codec.compress_file(
            in_path, out_path, ref_path=args.ref, threads=args.threads, try_xz=args.max
        )
        return "CRAM (reference-based)"
    elif fmt == "sam":
        sam_codec.compress_file(in_path, out_path, threads=args.threads, try_xz=args.max)
        return "columnar SAM (no reference given -- pass --ref for a bigger win)"
    else:
        generic_codec.compress_file(in_path, out_path, threads=args.threads, try_xz=args.max)
        return "generic"


def _decompress_one(in_path: Path, out_path: Path, args) -> None:
    format_id, _ = container.read_container(in_path)
    if format_id == container.FORMAT_FASTQ:
        fastq_codec.decompress_file(in_path, out_path)
    elif format_id == container.FORMAT_SAM:
        sam_codec.decompress_file(in_path, out_path, threads=args.threads)
    elif format_id == container.FORMAT_CRAM:
        cram_codec.decompress_file(in_path, out_path, ref_path=args.ref, threads=args.threads)
    else:
        generic_codec.decompress_file(in_path, out_path)


def cmd_compress(args):
    in_path = Path(args.input)

    if in_path.is_dir():
        out_dir = Path(args.output) if args.output else Path(str(in_path) + "_bioz")
        out_dir.mkdir(parents=True, exist_ok=True)
        files = sorted(p for p in in_path.iterdir() if p.is_file())
        if not files:
            print(f"no files found directly in {in_path}")
            return
        total_orig = total_comp = 0
        for f in files:
            out_path = out_dir / (f.name + ".bioz")
            t0 = time.time()
            try:
                mode = _compress_one(f, out_path, args)
            except Exception as e:
                print(f"  [skip] {f.name}: {e}")
                continue
            dt = time.time() - t0
            orig, comp = f.stat().st_size, out_path.stat().st_size
            ratio = orig / comp if comp else float("inf")
            total_orig += orig
            total_comp += comp
            print(f"  [{mode}] {f.name}: {_human(orig)} -> {_human(comp)} ({ratio:.2f}x, {dt:.1f}s)")
        overall = total_orig / total_comp if total_comp else float("inf")
        print(f"folder total: {_human(total_orig)} -> {_human(total_comp)} ({overall:.2f}x) -> {out_dir}/")
        return

    out_path = Path(args.output) if args.output else in_path.with_suffix(in_path.suffix + ".bioz")
    t0 = time.time()
    mode = _compress_one(in_path, out_path, args)
    dt = time.time() - t0
    real_orig = in_path.stat().st_size
    real_comp = out_path.stat().st_size
    ratio = real_orig / real_comp if real_comp else float("inf")
    print(f"[{mode}] {in_path.name}: {_human(real_orig)} -> {_human(real_comp)} "
          f"({ratio:.2f}x, {dt:.1f}s) -> {out_path}")


def cmd_decompress(args):
    in_path = Path(args.input)

    if in_path.is_dir():
        out_dir = Path(args.output) if args.output else Path(str(in_path) + "_decompressed")
        out_dir.mkdir(parents=True, exist_ok=True)
        files = sorted(p for p in in_path.iterdir() if p.is_file() and p.suffix == ".bioz")
        if not files:
            print(f"no .bioz files found directly in {in_path}")
            return
        for f in files:
            out_path = out_dir / f.stem  # strips the trailing .bioz
            t0 = time.time()
            try:
                _decompress_one(f, out_path, args)
            except Exception as e:
                print(f"  [skip] {f.name}: {e}")
                continue
            dt = time.time() - t0
            print(f"  {f.name} -> {out_path.name} ({dt:.1f}s)")
        print(f"folder total -> {out_dir}/")
        return

    if not args.output:
        out_path = in_path.with_suffix("") if in_path.suffix == ".bioz" else Path(str(in_path) + ".out")
    else:
        out_path = Path(args.output)

    t0 = time.time()
    _decompress_one(in_path, out_path, args)
    dt = time.time() - t0
    print(f"decompressed -> {out_path} ({dt:.1f}s)")


def main(argv=None):
    p = argparse.ArgumentParser(prog="bioz", description="Format-aware compressor for sequencing data (and anything else).")
    sub = p.add_subparsers(dest="cmd", required=True)

    pc = sub.add_parser("compress", help="compress a file, or every file in a folder, to .bioz")
    pc.add_argument("input", help="a file, or a folder -- each file inside is compressed with the codec matching its own type (fastq-aware, CRAM/columnar SAM, or generic)")
    pc.add_argument("-o", "--output", help="output path (file mode: default <input>.bioz; folder mode: default <input>_bioz/)")
    pc.add_argument("--generic", action="store_true", help="force generic mode for every file, skip fastq/sam-aware parsing")
    pc.add_argument("--lossy-quality", action="store_true", help="bin FASTQ quality scores (8 levels) for extra savings; sequence/ids/structure stay exactly lossless. OFF by default (fully lossless).")
    pc.add_argument("--ref", help="reference FASTA for BAM/SAM/CRAM input -> enables CRAM (reference-based) compression for every BAM/SAM found, the biggest win for aligned data. Without it (the default), falls back to columnar SAM splitting -- no reference is used unless you pass this.")
    pc.add_argument("--max", action="store_true", default=True, help="also try xz and keep whichever backend wins (default: on)")
    pc.add_argument("--fast", dest="max", action="store_false", help="zstd only, skip the xz comparison (faster, usually slightly worse ratio)")
    pc.add_argument("-T", "--threads", type=int, default=0, help="threads for backend compressors (default: all cores)")
    pc.set_defaults(func=cmd_compress)

    pd = sub.add_parser("decompress", help="decompress a .bioz file, or every .bioz file in a folder")
    pd.add_argument("input", help="a .bioz file, or a folder of .bioz files -- each is decompressed with the codec recorded in its own container header")
    pd.add_argument("-o", "--output", help="output path (file mode: default strips .bioz; folder mode: default <input>_decompressed/)")
    pd.add_argument("--ref", help="reference FASTA to decode a CRAM-format .bioz file (defaults to the path recorded at compress time, if still present)")
    pd.add_argument("-T", "--threads", type=int, default=0, help="threads for samtools/backend decompression (default: all cores)")
    pd.set_defaults(func=cmd_decompress)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
