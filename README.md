# bioz

Format-aware compressor. For FASTQ, it separates the file into homogeneous
streams (read IDs, `+` lines, 2-bit-packed sequence, quality scores) and
compresses each with the best of zstd/xz, instead of compressing the raw
interleaved 4-line text the way gzip/generic tools do. Any other file type
falls back to a generic best-of-zstd/xz pass.

Sequence, read IDs, and file structure are **always exactly lossless**
regardless of settings (verified by `tests/test_roundtrip.py`, including a
mixed-case/N-heavy edge case). Quality scores are lossless by default; pass
`--lossy-quality` to bin them to 8 levels for substantially better ratios
(this is the standard "quality binning" approach used elsewhere in the
field, not something exotic).

For SAM/BAM (aligned reads), pass `--ref` when you have the reference FASTA
used for alignment: this converts to CRAM, which needs the reference to
avoid storing redundant sequence -- the only technique here that actually
beat generic compression on BAM (which is already internally compressed,
so gzip/zstd/xz gain ~0% on it; CRAM still gets ~1.9x). Without `--ref`,
SAM/BAM falls back to columnar splitting (no reference needed, still beats
generic compression, just by less). `MD`/`NM` tags are regenerated (not
stored) by CRAM on decode since they're derivable from reference+CIGAR+seq
-- expected CRAM behavior, not data loss (verified in tests).

Both defaults are the safe ones: quality scores are **lossless** unless
you pass `--lossy-quality`, and SAM/BAM never uses a reference unless you
explicitly pass `--ref` (no auto-detection of reference files lying
around).

## Usage

```
bioz compress reads.fastq                    # -> reads.fastq.bioz, lossless
bioz compress reads.fastq --lossy-quality     # smaller, quality scores binned
bioz compress aligned.bam --ref genome.fa     # -> CRAM-based, needs the same --ref to decompress
bioz compress aligned.sam                     # no ref -> columnar SAM, still a real win
bioz compress anything.bin                    # generic path, any file type
bioz decompress reads.fastq.bioz              # -> reads.fastq (exact, if lossless)
bioz decompress aligned.bam.bioz --ref genome.fa -o aligned.bam
```

### Folder mode

Point `compress`/`decompress` at a folder instead of a single file, and
bioz compresses every file inside it, picking the right codec per file
by its own type -- a FASTQ gets fastq-aware lossless compression, a BAM
gets CRAM (only if you pass `--ref`) or columnar SAM otherwise, anything
else gets the generic fallback. One flat pass, not recursive into
subfolders.

```
bioz compress mydata/                         # -> mydata_bioz/<name>.bioz for each file
bioz compress mydata/ --ref genome.fa -o out/  # BAMs in the folder use CRAM against genome.fa
bioz decompress mydata_bioz/                   # -> mydata_bioz_decompressed/<name>, each via its own recorded codec
```

`--fast` skips the xz comparison (zstd only) for speed on large files.
`-T N` controls compressor thread count (default: all cores).

pod5 (raw nanopore signal) is not handled by bioz.

## Status / known limitation

Every codec (FASTQ, columnar SAM, CRAM, generic) streams through temp
files in bounded-size chunks rather than loading the whole input into
Python -- verified on real multi-GB files with `/usr/bin/time -v`: peak
memory stays roughly constant as input size grows (measured ~35MB for the
Python process itself regardless of file size, on files from a few MB up
to 1.7GB; see below). No fixed file-size ceiling from bioz's own side.

The one caveat: the *total* process-tree peak RSS reported by `time -v`
is higher than that (~1-2.7GB on the files tested) and grows somewhat
with input size -- that's `zstd --ultra -22 --long=27` (the backend's
strongest setting) running multi-threaded across all cores, an external
tool characteristic unrelated to bioz's own memory handling, not new in
this pass. If footprint at extreme scale (100GB+) matters more than
ratio, lower `-T`/pass `--fast` to trade some compression ratio for a
smaller, more predictable memory profile.

Everything else (correctness, ratio vs. gzip/zstd/xz) is tested and
benchmarked in `tests/`.

## Tests

```
python3 tests/make_test_data.py   # regenerate synthetic + real-data-slice test files
python3 tests/test_roundtrip.py   # correctness (must print ALL TESTS PASSED)
python3 tests/benchmark.py        # ratio/speed vs gzip -9 / zstd -19 / xz -9e
```
