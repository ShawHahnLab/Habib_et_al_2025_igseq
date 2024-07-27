#!/usr/bin/env python

"""
Aggregate information from a D gene initial DB, IgDiscover, and MINING-D.

This takes a D reference FASTA, an expressed_D.tab from IgDiscover, and an
output FASTA from MINING-D, and creates a CSV table summarizing evidence in
support of each D.  Every sequence in the reference FASTA gets at least one row
in the output, with individual rows for each detected match between the
reference and MINING-D and with the gene expression details repeated across
rows.  (Reference sequences with no match via either MINING-D or the expression
table will have blanks for those columns.)
"""

import argparse
from csv import DictReader, DictWriter
from collections import defaultdict
from Bio import SeqIO

FIELDS = [
    "gene", "gene_sequence", "candidate_sequence", "candidate_subsequence",
    "candidate_subsequence_fract","gene_subsequence_fract", "candidate_matches_total",
    "exp_count", "exp_unique_CDR3", "exp_unique_V", "exp_unique_J"]

def longest_submatch(query, ref):
    """Return the longest-matching subsequence of query in ref."""
    for sublen in range(len(query)+1)[::-1]:
        for idx_start in range(len(query)-sublen):
            idx_end = idx_start+sublen+1
            sub = query[idx_start:idx_end]
            if sub in ref:
                return sub
    return ""

def match_mining_d_to_ref(mining_d, d_ref):
    """Make dictionary matching MINING-D seq to reference seqs
    
    For each MINING-D sequence, this creates an inner dictionary mapping
    reference IDs to the matching MINING-D subsequences.  Only the longest
    subsequences are kept for any MINING-D sequence, but ties will result in
    multiple references noted per MINING-D sequence.
    """
    matches = defaultdict(dict)
    for seqid, seq in mining_d.items():
        for refid, ref in d_ref.items():
            sub = longest_submatch(seq, ref)
            matches[seqid][refid] = sub
    # reorder existing dicts and drop shorter matches
    for seqid, refs in matches.items():
        pairs = list(refs.items())
        pairs.sort(key = lambda row: -len(row[1]))
        maxlen = len(pairs[0][1])
        for pair in pairs:
            val = refs.pop(pair[0])
            if len(val) == maxlen:
                refs[pair[0]] = val
    return matches

def tabulate(d_ref, d_exp, mining_d, matches):
    """Combine parsed D seq details into consistent list of dictionaries"""
    # ref ID -> list of (D candidate, subsequence)
    matches_rev = defaultdict(set)
    for seqid, refs in matches.items():
        for refid, sub in refs.items():
            matches_rev[refid].add((seqid, sub))
    rows_out = []
    def handle_cands(cands):
        """Append rows for any MINING-D candidate matches"""
        for seqid, sub in cands:
            cand_matches_total = len(matches[seqid])
            cand_fract = len(sub)/len(mining_d[seqid])
            ref_fract = len(sub)/len(ref)
            row_out = row_out_template.copy()
            # No need to include the IDs from MINING-D, actually, since they're
            # arbitrary and there will be a 1:1 match with sequence content
            # anyway
            #row_out["candidate"] = seqid
            row_out["candidate_sequence"] = mining_d[seqid]
            row_out["candidate_subsequence"] = sub
            row_out["candidate_subsequence_fract"] = f"{cand_fract:.2f}"
            row_out["gene_subsequence_fract"] = f"{ref_fract:.2f}"
            row_out["candidate_matches_total"] = cand_matches_total
            rows_out.append(row_out)
    for refid, ref in d_ref.items():
        expression = d_exp.get(refid, {})
        row_out_template = {key: "" for key in FIELDS}
        row_out_template["gene"] = refid
        row_out_template["gene_sequence"] = ref
        row_out_template["exp_count"] = expression.get("count")
        row_out_template["exp_unique_CDR3"] = expression.get("unique_CDR3")
        row_out_template["exp_unique_V"] = expression.get("unique_V")
        row_out_template["exp_unique_J"] = expression.get("unique_J")
        cands = matches_rev[refid]
        if cands:
            handle_cands(cands)
        else:
            rows_out.append(row_out_template)
    rows_out.sort(key=lambda row: (row["gene"], row["candidate_matches_total"]))
    return rows_out

def aggregate_d(path_d_ref, path_d_exp, path_mining_d, path_out):
    """Aggregate information from a D gene initial DB, IgDiscover, and MINING-D"""
    # details of observed D expression
    d_exp = {}
    with open(path_d_exp, encoding="ASCII") as f_in:
        for row in DictReader(f_in, delimiter="\t"):
            d_exp[row["gene"]] = row
    # complete KIMDB D ref
    d_ref = {}
    for record in SeqIO.parse(path_d_ref, "fasta"):
        d_ref[record.id] = str(record.seq)
    # MINING-D D candidates
    mining_d = {}
    for record in SeqIO.parse(path_mining_d, "fasta"):
        mining_d[record.id] = str(record.seq)
    matches = match_mining_d_to_ref(mining_d, d_ref)
    rows_out = tabulate(d_ref, d_exp, mining_d, matches)
    with open(path_out, "w", encoding="ASCII") as f_out:
        writer = DictWriter(f_out, FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)

def main():
    """CLI for aggregate_d"""
    parser = argparse.ArgumentParser()
    parser.add_argument("d_ref", help="path to D reference FASTA")
    parser.add_argument("d_exp", help="path to IgDiscover output expressed_D.tab")
    parser.add_argument("mining_d", help="path to MINING-D output FASTA")
    parser.add_argument("output", help="path to output CSV")
    args = parser.parse_args()
    aggregate_d(args.d_ref, args.d_exp, args.mining_d, args.output)

if __name__ == "__main__":
    main()
