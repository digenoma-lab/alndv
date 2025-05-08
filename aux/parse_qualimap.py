#!/usr/bin/env python3
"""
parse_qualimap.py
-----------------
Extract key metrics from Qualimap BamQC genome_result(s).

Metrics collected
-----------------
* number_of_reads
* number_of_mapped_reads
* number_of_duplicated_reads
* median_insert_size
* mean_mapping_quality
* general_error_rate
* number_of_mismatches
* mean_coverage

Usage
-----
    # parse all reports in current directory
    ./parse_qualimap.py

    # or pass specific files / wildcards
    ./parse_qualimap.py sample*/genome_result*.txt
"""

import re, sys, glob, pathlib
from collections import defaultdict

# --------------------------------------------------------------------------- #
# Regex → (metric_key, capturing group index)
# --------------------------------------------------------------------------- #
PATTERNS = {
    "number_of_reads"         : re.compile(r"number of reads\s*=\s*([\d,]+)"),
    "number_of_mapped_reads"  : re.compile(r"number of mapped reads\s*=\s*([\d,]+)"),
    "number_of_duplicated_reads": re.compile(r"number of duplicated reads.*=\s*([\d,]+)"),
    "median_insert_size"      : re.compile(r"median insert size\s*=\s*([\d\.\,]+)"),
    "mean_mapping_quality"    : re.compile(r"mean mapping quality\s*=\s*([\d\.\,]+)"),
    "general_error_rate"      : re.compile(r"general error rate\s*=\s*([\d\.eE-]+)"),
    "number_of_mismatches"    : re.compile(r"number of mismatches\s*=\s*([\d,]+)"),
    "mean_coverage"           : re.compile(r"mean coverageData\s*=\s*([\d\.]+)"),
}

FIELDS = list(PATTERNS.keys())

def parse_report(path):
    """Return dict {metric: value} for one genome_result file."""
    txt = pathlib.Path(path).read_text(encoding="utf-8", errors="ignore")
    out = {}
    for key, pat in PATTERNS.items():
        m = pat.search(txt)
        if m:
            val = m.group(1).replace(",", "")
            # remove trailing X on mean_coverage if present
            if key == "mean_coverage":
                val = val.rstrip("Xx")
            out[key] = float(val) if "." in val or "e" in val.lower() else int(val)
        else:
            out[key] = "NA"
    return out

def main(files):
    # header
    print("file\t" + "\t".join(FIELDS))
    for fp in files:
        stats = parse_report(fp)
        values = [str(stats[f]) for f in FIELDS]
        #print(f"{pathlib.Path(fp).name}\t" + "\t".join(values))
        print(fp+"\t" + "\t".join(values))

if __name__ == "__main__":
    reports = sys.argv[1:] or glob.glob("*.genome_result*.txt") or glob.glob("*.genome_results*.txt")
    if not reports:
        sys.exit("No Qualimap genome_result.txt files found.")
    main(sorted(reports))

