#!/usr/bin/env python

"""
Make a FASTA for a figure S3 panel based on junction NT seqs

The logic of putting the sequences into alignments here is is *very* brittle,
but I'm just including it as a sanity-check on the example intermediate lineage
members shown in Figure S3.  See also: the "island" selection in the SONAR
module 2 snakemake rules, to manually select full sets of lineage members as we
actually did originally.
"""

import re
import sys
import csv
import argparse
from Bio import SeqIO

def _gap(seq, gap_pos, gap_len):
    return seq[:(gap_pos-1)] + "-"*gap_len + seq[(gap_pos-1):]

def _infer_subject(lineage):
    with open("metadata/isolates.csv", encoding="UTF8") as f_in:
        for row in csv.DictReader(f_in):
            if row["antibody_lineage"] == lineage:
                return row["subject"]
    raise ValueError

def _infer_specimen(subject):
    specimen = set()
    with open("metadata/biosamples.tsv", encoding="UTF8") as f_in:
        for row in csv.DictReader(f_in, delimiter="\t"):
            if row["igseq_Specimen_Subject"] == subject and \
                    row["igseq_Specimen_CellType"] == "IgG+":
                specimen.add(row["igseq_Specimen"])
    if len(specimen) != 1:
        raise ValueError
    return specimen.pop()

def _infer_isolate(lineage):
    with open("metadata/isolates.csv", encoding="UTF8") as f_in:
        for row in csv.DictReader(f_in):
            if row["antibody_lineage"] == lineage:
                return row
    raise ValueError

def _infer_junct_start(path, seq_id):
    for rec in SeqIO.parse(path, "fasta"):
        if rec.id == seq_id:
            return len(rec.seq) - 3
    raise ValueError(f"No cysTruncated {seq_id}?")

def _gather_members(path_junctions, path_seqs):
    """Get bulk NGS sequences with junctions matching the expected seqs"""
    # only keep representative examples as given in the this list, and use this
    # order for the output
    with open(path_junctions, encoding="ASCII") as f_in:
        juncts = [line.strip() for line in f_in]
    keepers = []
    juncts_left = juncts[:]
    for rec in SeqIO.parse(path_seqs, "fasta"):
        seq = str(rec.seq)
        junct = re.search("junction=([^ ]+)", rec.description).group(1)
        if junct in juncts_left:
            juncts_left.pop(juncts_left.index(junct))
            # truncate last NT
            keepers.append([rec.id, junct, seq[:-1]])
    keepers.sort(key=lambda info: juncts.index(info[1]))
    return keepers

def _make_germ_combo(germ):
    v_end   = re.search("-*$", germ["V"]).start() # first end gap
    d_start = re.match("-*", germ["D"]).end()     # first non-gap
    d_end   = re.search("-*$", germ["D"]).start() # first end gap
    j_start = re.match("-*", germ["J"]).end()     # first non-gap
    return germ["V"][:v_end] + \
        ("N" * max(0, d_start - v_end)) + \
        germ["D"][d_start:d_end] + \
        ("N" * max(0, j_start - d_end)) + \
        germ["J"][max(d_end, j_start):]

def _load_germline(path_germ, lineage):
    germline = {"V": {}, "D": {}, "J": {}}
    for segment, germ_seqs in germline.items():
        for rec in SeqIO.parse(f"{path_germ}/{segment}.fasta", "fasta"):
            germ_seqs[rec.id] = str(rec.seq)
    with open("analysis/isolates/summary_by_lineage.csv", encoding="UTF8") as f_in:
        for row in csv.DictReader(f_in):
            if row["antibody_lineage"] == lineage:
                break
        else:
            raise ValueError
    germ_here = {"V": row["vh"], "D": row["dh"], "J": row["jh"]}
    germ_seqs_here = {}
    for segment, germ_seqs in germline.items():
        germ_seqs_here[segment] = germ_seqs[germ_here[segment]]
    # truncate last NT
    germ_seqs_here["J"] = germ_seqs_here["J"][:-1]
    return germ_here, germ_seqs_here

def fig_s3(lineage, path_junctions, d_pos, path_out, gap_pos=0, gap_len=0):
    """Make a FASTA for a figure S3 panel based on junction NT seqs"""
    subject = _infer_subject(lineage)
    specimen = _infer_specimen(subject)
    isolate = _infer_isolate(lineage)
    germ_here, germ_seqs_here = _load_germline(f"analysis/germline-genbank/{subject}.IGH", lineage)
    junct_start = _infer_junct_start(
        f"analysis/sonar/{subject}.IGH/{specimen}/germline/V.cysTruncated.fasta",
        germ_here["V"])
    keepers = _gather_members(
        path_junctions,
        f"analysis/sonar/{subject}.IGH/{specimen}/output/sequences/"
            f"nucleotide/{specimen}_goodVJ_unique.fa")

    seqs = [info[2] for info in keepers] + [isolate["heavy_sequence"]]
    maxlen = max(len(s) for s in seqs)
    d_start = junct_start + d_pos - 1

    germ_seqs_here = {
        "V": germ_seqs_here["V"].ljust(maxlen, "-"),
        "D": _gap(germ_seqs_here["D"].rjust(
            d_start + len(germ_seqs_here["D"]), "-"), junct_start + gap_pos, gap_len).ljust(
                maxlen, "-"),
        "J": germ_seqs_here["J"].rjust(maxlen, "-")} 

    germ_combo = _make_germ_combo(germ_seqs_here)

    out = []
    for segment in ["V", "D", "J"]:
        out.append([germ_here[segment], germ_seqs_here[segment]])
    for info in keepers:
        if len(info[2].replace("-", "")) < maxlen:
            info[2] = _gap(info[2], junct_start + gap_pos, gap_len)
        out.append([info[0], info[2]])
    out.append([isolate["antibody_isolate"], isolate["heavy_sequence"]])
    out.append(["germline", germ_combo])

    with open(path_out, "w", encoding="ASCII") as f_out:
        for entry in out:
            f_out.write(">" + entry[0] + "\n" + entry[1] + "\n")

def main():
    """CLI for fig_s3"""
    parser = argparse.ArgumentParser()
    parser.add_argument("lineage")
    parser.add_argument("path_junctions")
    parser.add_argument("d_pos", type=int)
    parser.add_argument("path_out")
    parser.add_argument("--gap-pos", type=int, default=0)
    parser.add_argument("--gap-len", type=int, default=0)
    args = parser.parse_args()
    fig_s3(**vars(args))

if __name__ == "__main__":
    main()
