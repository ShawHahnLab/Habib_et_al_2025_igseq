#!/usr/bin/env python

"""
Convert GenBank flat file text into CSV of antibody isolate information.
"""

import re
import sys
import csv
import gzip
from collections import defaultdict
from Bio import SeqIO

def parse_gb_entry(gbf):
    for feat in gbf.features:
        if feat.type == "source":
            break
    else:
        raise ValueError(f"{gbf.name} no source?")
    match = re.match(
        r"Rhesus macaque (.*) antibody lineage (.*) antibody (.*) "
        r"isolated at ([0-9]+) weeks post-infection, (IG[HKL]) sequence",
        feat.qualifiers["note"][0])
    tissue = feat.qualifiers["tissue_type"][0]
    if tissue == "blood (PBMCs)":
        material = "PBMC"
    elif tissue == "lymph node tissue":
        material = "LN"
    else:
        raise ValueError(f"{gbf.name} {tissue}?")
    subject, lineage, isolate, timepoint, locus = match.groups()
    timepoint = int(timepoint)
    chain = "heavy" if locus == "IGH" else "light"
    seq_key = chain.capitalize() + "Seq"
    attrs_here = {
        "AntibodyIsolate": isolate,
        "AntibodyLineage": lineage,
        "Subject": subject,
        "Timepoint": timepoint,
        "Material": material,
        seq_key: gbf.seq
        }
    if chain == "light":
        attrs_here["LightLocus"] = locus
    return attrs_here

def convert_gb_isolates(gbs, csv_out):
    isolates = defaultdict(dict)
    for gb_path in gbs:
        with gzip.open(gb_path, "rt") as f_in:
            for gbf in SeqIO.parse(f_in, "gb"):
                attrs_here = parse_gb_entry(gbf)
                attrs = isolates[attrs_here["AntibodyIsolate"]]
                attrs.update(attrs_here)
    isolates = list(isolates.values())
    isolates.sort(
        key=lambda row: (row["AntibodyLineage"], row["AntibodyIsolate"], row["Timepoint"]))
    fields = list(isolates[0].keys())
    fields.sort(key=lambda txt: "Seq" in txt)
    with open(csv_out, "w") as f_out:
        writer = csv.DictWriter(f_out, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(isolates)

if __name__ == "__main__":
    convert_gb_isolates([sys.argv[1], sys.argv[2]], sys.argv[3])
