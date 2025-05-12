#!/usr/bin/env python3
"""
bedgz_to_matrix_autosomes_sex.py
  Merge multiple BED-4 gzip files into a coverage matrix
  **keeping only autosomes (1-22) plus X and Y**.
"""

import sys, gzip, csv
from collections import defaultdict

if len(sys.argv) < 2:
    sys.stderr.write(f"usage: {sys.argv[0]} file1.bed.gz [file2.bed.gz …]\n")
    sys.exit(1)

bed_paths  = sys.argv[1:]
sample_ids = [p.rsplit("/", 1)[-1].removesuffix(".depth.regions.bed.gz") for p in bed_paths]

# ------------------------------------------------------------------
# 1.  helpers
# ------------------------------------------------------------------
# valid chroms in both "1"/"chr1" style
_autosomes = {str(i) for i in range(1, 23)}
_sex       = {"X", "Y"}
VALID      = _autosomes | _sex

def canonical(chrom: str) -> str:
    """Return chromosome name without the leading 'chr' if present."""
    return chrom[3:] if chrom.lower().startswith("chr") else chrom

def keep(chrom: str) -> bool:
    """True if chrom is in the allowed set."""
    return canonical(chrom) in VALID

def sort_key(k):
    chrom, start, _ = k
    core = canonical(chrom)
    if core.isdigit():
        return (int(core), start)      # autosomes 1-22 numeric
    if core == "X":
        return (23, start)
    if core == "Y":
        return (24, start)
    return (99, start)                 # anything else goes to the end

# ------------------------------------------------------------------
# 2.  build matrix
# ------------------------------------------------------------------
matrix = defaultdict(lambda: [0.0] * len(bed_paths))

for idx, path in enumerate(bed_paths):
    with gzip.open(path, "rt") as fh:
        for chrom, start, end, cov in csv.reader(fh, delimiter="\t"):
            if not keep(chrom):
                continue               # skip unwanted chromosomes
            key = (chrom, int(start), int(end))
            matrix[key][idx] = float(cov)

rows = sorted(matrix.items(), key=lambda kv: sort_key(kv[0]))

# ------------------------------------------------------------------
# 3.  write TSV to stdout
# ------------------------------------------------------------------
out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
out.writerow(["chrom", "start", "end", *sample_ids])
for (chrom, start, end), covs in rows:
    out.writerow([chrom, start, end, *covs])

