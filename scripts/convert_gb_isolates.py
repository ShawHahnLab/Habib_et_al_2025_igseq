#!/usr/bin/env python

"""
Gather isolate CSVs into one unified table.
"""

import sys
import csv
from collections import defaultdict

def convert_gb_isolates(csvs_in, csv_out):
    isolates = defaultdict(dict)
    for csv_path in csvs_in:
        with open(csv_path) as f_in:
            for row in csv.DictReader(f_in):
                isolate = row["antibody_isolate"].replace("RHA1.V2", "RHA1")
                chain = row.get("chain", "heavy" if row.get("locus") == "IGH" else "light")
                seq_key = f"{chain}_sequence"
                attrs_here = {
                    "antibody_isolate": isolate,
                    "antibody_lineage": row["antibody_lineage"],
                    "subject": row["subject"],
                    "timepoint": int(row["timepoint"]),
                    seq_key: row["sequence"]}
                if chain == "light":
                    locus = "IGL" if row["subject"] == "5695" else row["locus"]
                    attrs_here["light_locus"] = locus
                attrs = isolates[attrs_here["antibody_isolate"]]
                attrs.update(attrs_here)
    isolates = list(isolates.values())
    isolates.sort(
        key=lambda row: (row["antibody_lineage"], row["antibody_isolate"], row["timepoint"]))
    fields = ["antibody_isolate", "antibody_lineage", "subject", "timepoint",
        "light_locus", "heavy_sequence", "light_sequence"]
    with open(csv_out, "w") as f_out:
        writer = csv.DictWriter(f_out, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(isolates)

if __name__ == "__main__":
    convert_gb_isolates([sys.argv[1], sys.argv[2]], sys.argv[3])
