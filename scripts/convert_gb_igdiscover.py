#!/usr/bin/env python

"""
Convert GenBank flat file text into CSV of IgDiscover information.
"""

import re
import sys
import csv
import gzip
from Bio import SeqIO

def parse_gb_entry(gbf):
    for feat in gbf.features:
        if feat.type == "source":
            break
    else:
        raise ValueError(f"{gbf.name} no source?")
    match = re.match(
        r"Rhesus Macaque (.*) antibody germline sequence for locus "
        r"(IG[HKL]) ([VJ]) segment inferred from antibody transcripts",
        feat.qualifiers["note"][0])
    subject, locus, segment = match.groups()
    seq_id = gbf.name.removeprefix(subject + "_")
    attrs = {
        "subject": subject,
        "sequence_id": seq_id,
        "segment": locus + segment,
        "sequence": gbf.seq}
    return attrs

def convert_gb_igdiscover(gb_path, csv_out):
    igdiscover = []
    with gzip.open(gb_path, "rt") as f_in:
        for gbf in SeqIO.parse(f_in, "gb"):
            igdiscover.append(parse_gb_entry(gbf))
    igdiscover.sort(key=lambda row: (row["subject"], row["segment"]))
    with open(csv_out, "w") as f_out:
        writer = csv.DictWriter(f_out, igdiscover[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(igdiscover)

if __name__ == "__main__":
    convert_gb_igdiscover(sys.argv[1], sys.argv[2])
