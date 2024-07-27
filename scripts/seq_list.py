#!/usr/bin/env python

"""
Convert FASTA to sorted list of sequences, one per line.

This is used here to get the same output from MINING-D for the same input after
disregarding sort order and the order-specific sequence IDs.
"""

import sys
from Bio import SeqIO

def seq_list(fasta_in, txt_out):
    with open(fasta_in) as f_in, open(txt_out, "w") as f_out:
        seqs = [str(record.seq) for record in SeqIO.parse(f_in, "fasta")]
        seqs.sort()
        for seq in seqs:
            f_out.write(f"{seq}\n")

if __name__ == "__main__":
    seq_list(sys.argv[1], sys.argv[2])
