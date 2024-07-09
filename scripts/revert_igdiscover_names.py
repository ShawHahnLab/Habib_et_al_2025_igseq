#!/usr/bin/env python

"""
IgDiscover sometimes makes a new name for a sequence even though it was already
in the starting database (when they are filtered out in one iteration and then
re-discovered in a later iteration, maybe?).  This switches IDs back when the
sequences match and the ID is derived from the ID in the starting database.

(This only handles the simple case where the sequence for <some ID>_S#### from
IgDiscover matches that of <some ID> in the starting database; there are other
instances where the sequence matches but names are different but I haven't
dealt with those.)
"""

import re
import sys
from Bio import SeqIO

def revert_igdiscover_names(ref_in, igdisc_in, igdisc_out):
    """Switch IgDiscover sequence IDs back for those matching a reference"""
    ref = {str(rec.seq): rec.id for rec in SeqIO.parse(ref_in, "fasta")}
    with open(igdisc_out, "w", encoding="ASCII") as f_out:
        for rec in SeqIO.parse(igdisc_in, "fasta"):
            orig_id = ref.get(str(rec.seq))
            if "_S" in rec.id and re.sub("_S.*", "", rec.id) == orig_id:
                sys.stderr.write(f"{rec.id} -> {orig_id}\n")
                rec.id = orig_id
            SeqIO.write(rec, f_out, "fasta-2line")

if __name__ == "__main__":
    revert_igdiscover_names(sys.argv[1], sys.argv[2], sys.argv[3])
