#!/usr/bin/env python

"""
Gather a set of FASTA files into a CSV.

The sequence will be included in a "sequence" column.  Use named capture groups
in the pattern to specify what fields map to which other columns:

    https://docs.python.org/3.9/library/re.html#index-17
"""

import re
import csv
import gzip
import argparse
from pathlib import Path
from Bio import SeqIO

OPENERS = {".gz": gzip.open}

def _parse_fields(txts):
    shared_vals = {}
    if txts:
        for pair in txts:
            key, val = pair.split("=", 1)
            shared_vals[key] = val
    return shared_vals

def _get_source_note(gbf):
    for feat in gbf.features:
        if feat.type == "source":
            break
    else:
        raise ValueError(f"{gbf.name} no source?")
    return feat.qualifiers["note"][0]

def tabulate_seqs(paths, output, pattern="", extra_field=None, fmt="fasta"):
    """
    paths: list of paths to FASTA files
    output: path to output CSV
    pattern: regular expression using Python named capture groups
    extra_fields: list of key=value pairs to use for constant columns
    fmt: BioPython input format ("fasta", "gb", ...)
    """
    shared_vals = _parse_fields(extra_field)
    with open(output, "w", encoding="UTF8") as f_out:
        writer = None
        for path in paths:
            path = Path(path)
            with OPENERS.get(path.suffix.lower(), open)(path, "rt", encoding="UTF8") as f_in:
                for rec in SeqIO.parse(f_in, fmt):
                    row = shared_vals.copy()
                    if pattern:
                        text = rec.description
                        if fmt == "gb":
                            text = _get_source_note(rec)
                        match = re.match(pattern, text)
                        if not match:
                            raise ValueError(
                                f'No match for "{text}" with pattern "{pattern}"')
                        row.update(match.groupdict())
                    if "sequence_id" not in row:
                        row["sequence_id"] = rec.id
                    row["sequence"] = str(rec.seq)
                    if not writer:
                        writer = csv.DictWriter(
                            f_out,
                            fieldnames=row.keys(),
                            lineterminator="\n")
                        writer.writeheader()
                    writer.writerow(row)

def main():
    """CLI for tabulate_seqs"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p", "--pattern", help="regex for parsing fields from seq descriptions")
    parser.add_argument(
        "-x", "--extra-field", action="append", help="key=value to apply for all rows")
    parser.add_argument(
        "-o", "--output", required="true", help="CSV output path")
    parser.add_argument(
        "-f", "--fmt", default="fasta", help="input format")
    parser.add_argument(
        "paths", nargs="*", help="input FASTA paths")
    args = parser.parse_args()
    tabulate_seqs(**vars(args))

if __name__ == "__main__":
    main()
