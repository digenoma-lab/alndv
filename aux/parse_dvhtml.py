#!/usr/bin/env python3
"""
extract_dv_counts.py
--------------------
Pull Biallelic_SNP, Biallelic_Deletion and Biallelic_Insertion counts
from one or many DeepVariant *.visual_report.html files.

USAGE
-----
  # scan every report in the current directory
  ./extract_dv_counts.py

  # or pass explicit files / wildcards
  ./extract_dv_counts.py sample1.html trio*.html

OUTPUT
------
TSV on stdout:

file	Biallelic_SNP	Biallelic_Deletion	Biallelic_Insertion
F50P1.autosomes.visual_report.html	3902365	382869	354257
...
"""

import re, sys, glob, pathlib
from html import unescape
from collections import defaultdict

# ────────────────────────────────────────────────────────────────────────────
# Regexes that cover both storage styles
# ────────────────────────────────────────────────────────────────────────────
# 1)  {"label": "Biallelic_SNP", "value": 3902365}
LBL_RE  = re.compile(r'"label"\s*:\s*"(?P<label>Biallelic_(?:SNP|Deletion|Insertion))"\s*,\s*"value"\s*:\s*(?P<val>\d+)', re.I)

# 2)  "Biallelic_SNP": 3902365
KV_RE   = re.compile(r'"(?P<label>Biallelic_(?:SNP|Deletion|Insertion))"\s*:\s*(?P<val>\d+)', re.I)

def parse_html(path):
    """Return dict with counts from one HTML report."""
    txt = pathlib.Path(path).read_text(encoding="utf-8", errors="ignore")
    txt = unescape(txt)             # unescape &quot; etc. just in case

    counts = defaultdict(int)

    # match style 1
    for m in LBL_RE.finditer(txt):
        counts[m['label']] = int(m['val'])

    # match style 2  (only if still missing)
    if len(counts) < 3:
        for m in KV_RE.finditer(txt):
            counts[m['label']] = int(m['val'])

    return counts

def main(files):
    print("file\tBiallelic_SNP\tBiallelic_Deletion\tBiallelic_Insertion")
    for fp in files:
        c = parse_html(fp)
        print(f"{pathlib.Path(fp).name}\t"
              f"{c.get('Biallelic_SNP',0)}\t"
              f"{c.get('Biallelic_Deletion',0)}\t"
              f"{c.get('Biallelic_Insertion',0)}")

if __name__ == "__main__":
    html_files = sys.argv[1:] or glob.glob("*.visual_report.html")
    if not html_files:
        sys.exit("No DeepVariant *.visual_report.html files found.")
    main(sorted(html_files))

