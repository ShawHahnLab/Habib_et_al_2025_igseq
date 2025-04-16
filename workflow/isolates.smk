### Isolates

def make_target_igblast_isolates():
    attrs = set()
    for row in METADATA["isolates"]:
        attrs.add((row["subject"], "IGH"))
        if row["light_sequence"]:
            attrs.add((row["subject"], row["light_locus"]))
    attrs = list(attrs)
    attrs.sort()
    attrs = {key: val for key, val in zip(("subject", "locus"), zip(*attrs))}
    return expand("analysis/isolates/{subject}.{locus}/igblast.tsv", zip, **attrs)

TARGET_IGBLAST_ISOLATES = make_target_igblast_isolates()

rule igblast_isolates_combined:
    """Combined TSV of IgBLAST AIRR across heavy+light seqs of all isolates"""
    output: "analysis/isolates/igblast.tsv"
    input: TARGET_IGBLAST_ISOLATES
    run:
        with open(output[0], "w") as f_out:
            writer = None
            for path in input:
                with open(path) as f_in:
                    reader = csv.DictReader(f_in, delimiter="\t")
                    for row in reader:
                        if not writer:
                            writer = csv.DictWriter(
                                f_out, reader.fieldnames, delimiter="\t", lineterminator="\n")
                            writer.writeheader()
                        writer.writerow(row)

rule igblast_isolates:
    """IgBLAST AIRR for individual heavy+light sequences for each isolate"""
    output: "analysis/isolates/{subject}.{locus}/igblast.tsv"
    input:
        query="analysis/isolates/{subject}.{locus}/query.fasta",
        ref=expand("analysis/{germline}/{{subject}}.{{locus}}/{segment}.fasta", germline=config["germline"], segment=["V", "D", "J"])
    params:
        ref=lambda w, input: Path(input.ref[0]).parent
    conda: "igseq.yml"
    shell: "igseq igblast -S rhesus -r {params.ref} -Q {input.query} -outfmt 19 -out {output}"

rule igblast_isolates_input:
    output: temp("analysis/isolates/{subject}.{locus}/query.fasta")
    input: "metadata/isolates.csv"
    run:
        with open(input[0]) as f_in, open(output[0], "w") as f_out:
            col_seq = "heavy_sequence" if wildcards.locus == "IGH" else "light_sequence"
            for row in csv.DictReader(f_in):
                if row["subject"] == wildcards.subject and \
                        wildcards.locus in (row["light_locus"], "IGH"):
                    seqid = row["antibody_isolate"]
                    seq = row[col_seq]
                    f_out.write(f">{seqid}\n{seq}\n")

rule igblast_isolates_lineage_summary:
    """Rough summary table of top gene calls by lineage according to IgBLAST"""
    output: "analysis/isolates/summary_by_lineage.csv"
    input: "analysis/isolates/igblast.tsv"
    run:
        from collections import defaultdict
        isolate_map = {row["antibody_isolate"]: row for row in METADATA["isolates"]}
        # gather up a tally by lineage by locus+segment of all observed IgBLAST
        # gene calls, counting entries from ties separately
        # (lineage -> locus+segment -> tally of gene calls across antibodies)
        calls = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        with open(input[0]) as f_in:
            for row in csv.DictReader(f_in, delimiter="\t"):
                lineage = isolate_map[row["sequence_id"]]["antibody_lineage"]
                for seg in ["v", "d", "j"]:
                    calls_here = row[f"{seg}_call"]
                    if calls_here and not "," in calls_here:
                        locseg = calls_here[:4]
                        for call in calls_here.split(","):
                            calls[lineage][locseg][call] += 1
        out = []
        colmap = {
            "IGHV": "vh", "IGHD": "dh", "IGHJ": "jh",
            "IGKV": "vl", "IGKJ": "jl",
            "IGLV": "vl", "IGLJ": "jl"}
        for lineage, locseg_calls in calls.items():
            tops = {}
            for locseg, calls_here in locseg_calls.items():
                # We went over all this lineage by lineage with a fine-toothed
                # comb manually so it's all a bit rough to just take IgBLAST's
                # output as-is here; with that in mind, just shunt D calls to
                # D3-15 as we know they should be, assuming that call at least
                # shows up in the list (which it better)
                if "IGHD3-15*01" in calls_here:
                    calls_here = "IGHD3-15*01"
                elif locseg == "IGHD":
                    raise ValueError("No D3-15?!")
                else:
                    calls_here = [(num, name) for name, num in calls_here.items()]
                    calls_here = sorted(calls_here)[0][1]
                tops[locseg] = calls_here
            row_out = {"antibody_lineage": lineage}
            for locseg, call in tops.items():
                colname = colmap[locseg]
                if colname in row_out:
                    # basically just to make sure no IGK/IGL collisions for any
                    # one lineage which would make no sense
                    raise ValueError(f"{row_out} but we also have {locseg}?")
                row_out[colname] = tops[locseg]
            out.append(row_out)
        with open(output[0], "w") as f_out:
            writer = csv.DictWriter(
                f_out,
                ["antibody_lineage", "vh", "dh", "jh", "vl", "jl"],
                lineterminator="\n")
            writer.writeheader()
            writer.writerows(out)

rule igblast_isolates_summary:
    """Per-isolate heavy+light summary table including some lineage+subject info"""
    output: "analysis/isolates/summary.csv"
    input:
        by_lineage="analysis/isolates/summary_by_lineage.csv",
        igblast="analysis/isolates/igblast.tsv"
    run:
        isolate_map = {row["antibody_isolate"]: row for row in METADATA["isolates"]}
        with open(input.by_lineage) as f_in:
            lineage_calls = {row["antibody_lineage"]: row for row in csv.DictReader(f_in)}
        out = defaultdict(dict)
        with open(input.igblast) as f_in:
            for row in csv.DictReader(f_in, delimiter="\t"):
                isolate_attrs = isolate_map[row["sequence_id"]]
                lineage = isolate_attrs["antibody_lineage"]
                subject = isolate_attrs["subject"]
                row_out = out[row["sequence_id"]]
                row_out["sequence_id"] = row["sequence_id"]
                row_out["antibody_lineage"] = lineage
                row_out["subject"] = subject
                row_out["timepoint"] = isolate_attrs["timepoint"]
                suffix = "h" if row["v_call"].startswith("IGH") else "l"
                for seg in ["v", "d", "j"]:
                    key = seg + suffix
                    if key in lineage_calls[lineage]:
                        row_out[key] = lineage_calls[lineage][key]
                row_out[f"shm_v{suffix}"] = round(100-float(row["v_identity"]), 1)
                if suffix == "h":
                    row_out["heavy_junction_aa"] = row["junction_aa"]
                    row_out["cdrh3_len"] = len(row["cdr3_aa"])
        out = list(out.values())
        with open(output[0], "w") as f_out:
            writer = csv.DictWriter(
                f_out,
                ["sequence_id", "antibody_lineage", "subject", "timepoint", "heavy_junction_aa",
                    "cdrh3_len", "vh", "dh", "jh", "vl", "jl", "shm_vh", "shm_vl"],
                lineterminator="\n")
            writer.writeheader()
            writer.writerows(out)
