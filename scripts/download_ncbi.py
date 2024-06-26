#!/usr/bin/env python

"""
Simple NCBI FASTA accession downloader.
"""

import os
import sys
import time
from io import StringIO
from urllib.error import HTTPError
from Bio import Entrez
from Bio import SeqIO
from subprocess import run

Entrez.email = os.getenv("ENTREZ_EMAIL")
if not Entrez.email:
    proc = run(["git", "config", "user.email"], capture_output=True, check=False, text=True)
    if proc.returncode:
        raise Exception("set ENTREZ_EMAIL environment variable to download from GenBank")
    Entrez.email = proc.stdout.strip()

def _repeating_efetch(**kwargs):
    while True:
        try:
            handle = Entrez.efetch(**kwargs)
        except HTTPError as error:
            sys.stderr.write(str(error))
            sys.stderr.write("\n")
            time.sleep(5)
        else:
            break
    return handle

def download_ncbi(db, acc, path_out=None, rettype="fasta"):
    """Download a FASTA file for a single accession.

    For GenBank, use db="nucleotide" and rettype "fasta" or "gb"
    For PDB, use db="protein"
    """
    handle = _repeating_efetch(db=db, id=acc, rettype=rettype, retmode="text")
    txt = StringIO(handle.read().strip() + "\n")
    fmt = "gb" if rettype == "gb" else "fasta-2line"
    recs = 0
    for record in SeqIO.parse(txt, rettype):
        if path_out:
            with open(path_out, "w") as f_out:
                SeqIO.write(record, f_out, fmt)
        else:
            SeqIO.write(record, sys.stdout, fmt)
        recs += 1
    if not recs:
        raise ValueError(f"No records found for {acc} (DB: {db}, rettype: {rettype})")

if __name__ == "__main__":
    download_ncbi(*sys.argv[1:])
