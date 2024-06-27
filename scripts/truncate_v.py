#!/usr/bin/env python

"""
Truncate V at the CDRH3 conserved cysteine for SONAR's ID/DIV calculations
"""

import sys
import csv
import subprocess
from io import StringIO

def truncate_v(fasta_in, fasta_out):
    """Truncate V at the conserved cysteine"""
    # (We can just let igseq default to the builtin rhesus macaque references
    # since that'll be plenty accurate enough to just report where FWR3 ends)
    result = subprocess.run(
        ["igseq", "igblast", "-S", "rhesus", "-Q", fasta_in, "-outfmt", "19"],
        stdout=subprocess.PIPE, text=True, check=True)
    with open(fasta_out, "w", encoding="ASCII") as f_out:
        for row in csv.DictReader(StringIO(result.stdout), delimiter="\t"):
            seqid = row["sequence_id"]
            pos = int(row["fwr3_end"])
            seq = row["sequence"][:pos]
            f_out.write(f">{seqid}\n{seq}\n")

if __name__ == "__main__":
    truncate_v(sys.argv[1], sys.argv[2])
