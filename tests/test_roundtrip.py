"""Correctness tests: compress -> decompress must reproduce the original
exactly (lossless mode), or reproduce everything except quality scores
exactly (lossy-quality mode). Run from anywhere: python3 tests/test_roundtrip.py
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
DATA = HERE / "test_data"
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"cmd failed: {cmd}\nstdout={r.stdout}\nstderr={r.stderr}")
    return r


def read_fastq_records(path):
    lines = Path(path).read_text().splitlines()
    return [tuple(lines[i : i + 4]) for i in range(0, len(lines), 4)]


failures = []


def check(name, cond):
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {name}")
    if not cond:
        failures.append(name)


# --- lossless fastq round-trip: must be byte-identical ---
for fq in ["small.fastq", "medium_shortreads.fastq", "real_slice.fastq"]:
    src = DATA / fq
    archive = WORK / (fq + ".bioz")
    out = WORK / (fq + ".roundtrip")
    run(["bioz", "compress", str(src), "-o", str(archive)])
    run(["bioz", "decompress", str(archive), "-o", str(out)])
    check(f"lossless round-trip byte-identical: {fq}", src.read_bytes() == out.read_bytes())

# --- lossy-quality fastq round-trip: ids/seq/plus identical, quality may differ ---
fq = "medium_shortreads.fastq"
src = DATA / fq
archive = WORK / (fq + ".lossy.bioz")
out = WORK / (fq + ".lossy.roundtrip")
run(["bioz", "compress", str(src), "-o", str(archive), "--lossy-quality"])
run(["bioz", "decompress", str(archive), "-o", str(out)])
orig_records = read_fastq_records(src)
new_records = read_fastq_records(out)
check("lossy: same number of records", len(orig_records) == len(new_records))
ids_match = all(o[0] == n[0] and o[1] == n[1] and o[2] == n[2] for o, n in zip(orig_records, new_records))
check("lossy: ids/sequence/plus lines exactly preserved", ids_match)
lens_match = all(len(o[3]) == len(n[3]) for o, n in zip(orig_records, new_records))
check("lossy: quality line lengths preserved", lens_match)
qual_differs_somewhere = any(o[3] != n[3] for o, n in zip(orig_records, new_records))
check("lossy: quality scores actually got binned (not silently lossless)", qual_differs_somewhere)
distinct_qual_bytes = len(set("".join(n[3] for n in new_records)))
check(f"lossy: quality alphabet collapsed to <=8 symbols (got {distinct_qual_bytes})", distinct_qual_bytes <= 8)

# --- generic (non-fastq) file round-trip ---
gfile = DATA / "generic.txt"
archive = WORK / "generic.txt.bioz"
out = WORK / "generic.txt.roundtrip"
run(["bioz", "compress", str(gfile), "-o", str(archive)])
run(["bioz", "decompress", str(archive), "-o", str(out)])
check("generic file round-trip byte-identical", gfile.read_bytes() == out.read_bytes())

# --- edge case: empty-ish / N-heavy sequences handled without crashing ---
edge = WORK / "edge.fastq"
edge.write_text(
    "@r1\n"
    "NNNNNNNNNN\n"
    "+\n"
    "!!!!!!!!!!\n"
    "@r2\n"
    "ACGTNacgtN\n"
    "+\n"
    "IIIIIIIIII\n"
)
archive = WORK / "edge.fastq.bioz"
out = WORK / "edge.fastq.roundtrip"
run(["bioz", "compress", str(edge), "-o", str(archive)])
run(["bioz", "decompress", str(archive), "-o", str(out)])
check("edge case (N-heavy / mixed-case bases) round-trip exact", edge.read_bytes() == out.read_bytes())

# --- SAM (columnar, no reference): .sam input should round-trip byte-exact ---
sam_src = DATA / "slice.sam"
sam_archive = WORK / "slice.sam.bioz"
sam_out = WORK / "slice.roundtrip.sam"  # must end in .sam so sam_codec knows to emit text, not BAM
run(["bioz", "compress", str(sam_src), "-o", str(sam_archive)])
run(["bioz", "decompress", str(sam_archive), "-o", str(sam_out)])
check("SAM (no ref) round-trip byte-identical", sam_src.read_bytes() == sam_out.read_bytes())

# --- SAM (columnar, no reference): .bam input -> record-level equivalence ---
bam_src = DATA / "slice.bam"
bam_archive = WORK / "slice.bam.nosam.bioz"
bam_out_sam = WORK / "slice.bam.nosam.roundtrip.sam"
run(["bioz", "compress", str(bam_src), "-o", str(bam_archive)])
run(["bioz", "decompress", str(bam_archive), "-o", str(bam_out_sam)])
expected_sam = run(["samtools", "view", "-h", str(bam_src)]).stdout
check("BAM (no ref) -> SAM record-level identical to samtools view", expected_sam.encode() == bam_out_sam.read_bytes())

# --- CRAM (reference-based): .bam input -> decompress to .bam -> record-level equivalence ---
ref = Path("/data/Raphael/DNAscentTest_I/sacCer3.fa")
if ref.exists():
    cram_archive = WORK / "slice.bam.cram.bioz"
    cram_out_bam = WORK / "slice.cram.roundtrip.bam"
    run(["bioz", "compress", str(bam_src), "-o", str(cram_archive), "--ref", str(ref)])
    run(["bioz", "decompress", str(cram_archive), "-o", str(cram_out_bam), "--ref", str(ref)])
    expected = run(["samtools", "view", str(bam_src)]).stdout.splitlines()
    got = run(["samtools", "view", str(cram_out_bam)]).stdout.splitlines()
    # Core alignment fields (QNAME..QUAL, the first 11 columns) must match exactly.
    # Optional tags are NOT compared byte-for-byte: CRAM deliberately doesn't store
    # MD/NM (they're recomputable from the reference + CIGAR + sequence) and
    # regenerates them on decode, which can reorder/add tags -- expected CRAM
    # behavior, not data loss. NM's *value* is checked instead, since that's the
    # one optional tag actually worth verifying semantically.
    core_match = len(expected) == len(got) and all(
        e.split("\t")[:11] == g.split("\t")[:11] for e, g in zip(expected, got)
    )
    check("CRAM round-trip: core alignment fields (QNAME..QUAL) identical", core_match)

    def nm_value(line):
        for f in line.split("\t")[11:]:
            if f.startswith("NM:i:"):
                return f
        return None

    nm_match = all(nm_value(e) == nm_value(g) for e, g in zip(expected, got))
    check("CRAM round-trip: NM tag values preserved", nm_match)

    # wrong-reference safety check: must refuse, not silently misdecode
    wrong_ref = WORK / "wrong_ref.fa"
    wrong_ref.write_text(">fake\nACGTACGTACGT\n")
    r = subprocess.run(["bioz", "decompress", str(cram_archive), "-o", str(WORK / "should_fail.bam"), "--ref", str(wrong_ref)],
                        capture_output=True, text=True)
    check("CRAM refuses to decode against a mismatched reference", r.returncode != 0)
else:
    print(f"[SKIP] CRAM tests (reference not found at {ref})")

print()
if failures:
    print(f"{len(failures)} FAILURE(S):", failures)
    sys.exit(1)
else:
    print("ALL TESTS PASSED")
