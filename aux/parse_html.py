

#!/usr/bin/env python3
"""
extract_dv_variant_counts.py
----------------------------
Parse DeepVariant visual‑report HTML files and grab the counts for
  • Biallelic_SNP
  • Biallelic_Deletion
  • Biallelic_Insertion

Usage
-----
    # all *.visual_report.html in the cwd
    python extract_dv_variant_counts.py

    # or specify files / wildcards explicitly
    python extract_dv_variant_counts.py F50P1*.html sample2.html

Output
------
TSV to stdout:

    file	Biallelic_SNP	Biallelic_Deletion	Biallelic_Insertion
    F50P1.autosomes.visual_report.html	3902365	382869	354257
    ...

Redirect to a file if you prefer:
    python extract_dv_variant_counts.py > dv_variant_counts.tsv
"""
import sys, re, glob, pathlib, json, html
from collections import defaultdict

# --------------------------------------------------------------------------- #
# Helper: harvest the JSON blob DeepVariant embeds in a <script> tag
# --------------------------------------------------------------------------- #
VAR_RE = re.compile(r'"Biallelic_(SNP|Deletion|Insertion)"\s*:\s*(\d+)', re.I)

def parse_file(path):
    counts = defaultdict(int)
    with open(path, encoding="utf‑8") as fh:
        txt = fh.read()
    # 1) fast regex grab (works for v1.4 + reports)
    for kind, num in VAR_RE.findall(txt):
        counts[kind.capitalize()] = int(num)
    # 2) fallback: look for the JSON inside  id="variant-type-data"
    if not counts:
        m = re.search(r'<script[^>]+id="variant-type-data"[^>]*>(.*?)</script>',
                      txt, re.S|re.I)
        if m:
            data = json.loads(html.unescape(m.group(1)))
            for k in ("Biallelic_SNP", "Biallelic_Deletion", "Biallelic_Insertion"):
                counts[k.split('_')[1].capitalize()] = int(data.get(k, 0))
    return counts

# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    files = sys.argv[1:] or glob.glob("*.visual_report.html")
    if not files:
        sys.exit("No HTML files found ‑ please give filenames or run in the report directory.")
    # header
    print("file\tBiallelic_SNP\tBiallelic_Deletion\tBiallelic_Insertion")
    for fp in files:
        counts = parse_file(fp)
        snp  = counts.get("Snp",        0)
        del_ = counts.get("Deletion",   0)
        ins  = counts.get("Insertion",  0)
        print(f"{pathlib.Path(fp).name}\t{snp}\t{del_}\t{ins}")

