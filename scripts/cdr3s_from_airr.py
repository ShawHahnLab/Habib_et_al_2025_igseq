#!/usr/bin/env python

import sys
import gzip
from collections import defaultdict

def cdr3s_from_airr(tsvgzs_in, fasta_out):
    seqs = defaultdict(int)
    for path in tsvgzs_in:
        opener = open
        if str(path).endswith(".gz"):
            opener = gzip.open
        with opener(path, "rt") as f_in:
            reader = DictReader(f_in, delimiter="\t")
            for row in reader:
                # If this is SONAR AIRR only use the cluster centroids
                if "cluster_count" in row and not row["cluster_count"]:
                    continue
                # IgBLAST puts both cdr3 and junction cols in its AIRR TSV
                # output while SONAR only provides junction.  We'll work
                # with either.
                if row["locus"] == "IGH" and (row.get("cdr3") or row.get("junction")):
                    try:
                        seq = row["cdr3"]
                    except KeyError:
                        seq = row["junction"][3:-3]
                    seqs[seq] += 1
    # Write deduplicated sequences sorted by abundance
    seqs = list(seqs.items())
    seqs.sort(key=lambda pair: -pair[1])
    with open(fasta_out, "wt") as f_out:
        for idx, pair in enumerate(seqs):
            seqid = f"cdr3_{idx}"
            desc = f"count={pair[1]}"
            seq = pair[0]
            f_out.write(f">{seqid} {desc}\n{seq}\n")

def main():
    cdr3s_from_airr(sys.argv[1:-1], sys.argv[-1])
