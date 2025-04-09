### Output
#
# Some final summary outputs approximating what's shown in the structure paper itself.

def final_gene_name(txt):
    """Format the gene names like the paper has them"""
    return re.sub(r"IG([HKL])([VDJ])", r"\2\1", txt)

def indexish(items, item, default=-1):
    """Index of item in items or -1 if not present (just to help override sorting)"""
    try:
        return items.index(item)
    except ValueError:
        return default

rule structure_output_fig1b:
    """Summary table of antibody sequence characteristics"""
    output: "output/structure_fig1b.csv"
    input: "analysis/isolates/summary.csv"
    run:
        mabs = [
            "6070-a.01", "42056-a.01", "5695-b.01", "T646-a.01", "41328-a.01", "V033-a.01",
            "44715-a.01", "40591-a.01", "6561-a.01", "42056-b.01", "V031-a.01"]
        out = []
        with open(input[0]) as f_in:
            # considered trying to be cute and automatically note the indel
            # status via the IgBLAST results too, but then remembered IgBLAST
            # is pretty terrible at figuring that out, especially midway
            # through a lineage like these are.  So, skipping that column.
            for row in csv.DictReader(f_in):
                if row["sequence_id"] in mabs:
                    out.append({
                        "mAb ID": row["sequence_id"],
                        "Macmul VH gene": final_gene_name(row["vh"]),
                        "Macmul VL gene": final_gene_name(row["vl"]),
                        "SHM IGHV": final_gene_name(row["shm_vh"]),
                        "SHM IGLV": final_gene_name(row["shm_vl"]),
                        "HCDR3 length (aa)": row["cdrh3_len"],
                        "VDJ Junction": row["heavy_junction_aa"]})
        out.sort(key=lambda row: mabs.index(row["mAb ID"]))
        with open(output[0], "w") as f_out:
            writer = csv.DictWriter(f_out, out[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(out)

rule structure_output_fig2a:
    """Summary table of antibody heavy chain germline sequence usage"""
    output: "output/structure_fig2a.csv"
    input: "analysis/isolates/summary.csv"
    run:
        mabs = [
            "42056-b.01", "6561-a.01", "40591-a.01", "T646-a.01", "V031-a.01",
            "6070-a.01", "5695-b.01", "RHA1.01", "44715-a.01", "41328-a.01",
            "42056-a.01", "V033-a.01"]
        out = []
        with open(input[0]) as f_in, open(output[0], "w") as f_out:
            writer = csv.DictWriter(
                f_out, ["Antibody ID", "Macmul VH gene", "Macmul DH gene", "Macmul JH gene"],
                lineterminator="\n")
            writer.writeheader()
            for row in csv.DictReader(f_in):
                if row["sequence_id"] in mabs:
                    out.append({
                        "Antibody ID": row["sequence_id"],
                        "Macmul VH gene": final_gene_name(row["vh"]),
                        "Macmul DH gene": final_gene_name(row["dh"]),
                        "Macmul JH gene": final_gene_name(row["jh"])})
            out.sort(key=lambda row: mabs.index(row["Antibody ID"]))
            writer.writerows(out)

rule structure_output_tableS2:
    """Summary table of animal, antibody, and lineage information"""
    output: "output/structure_tableS2.csv"
    input: "analysis/isolates/summary.csv"
    run:
        lineages = ["6070-a", "42056-a", "42056-b", "5695-b", "T646-a",
            "41328-a", "V033-a", "44715-a", "40591-a", "6561-a", "V031-a"]
        isolates = ["5695-b.04", "5695-b.05", "5695-b.02", "5695-b.03", "5695-b.01"]
        out = []
        with open(input[0]) as f_in:
            for row in csv.DictReader(f_in):
                splat = lambda keys: " ".join(final_gene_name(row[k]) for k in keys if row[k])
                if row["antibody_lineage"] in lineages:
                    out.append({
                        "Animal ID": "RM" + row["subject"],
                        "Timepoint": "wk" + row["timepoint"],
                        "mAb ID": row["sequence_id"],
                        "Heavy Ig genes": splat(("vh", "dh", "jh")),
                        "Light Ig genes": splat(("vl", "jl")),
                        "VH %nt mut": row["shm_vh"],
                        "VL %nt mut": row["shm_vl"],
                        "_lineage": row["antibody_lineage"]})
        def sorter(row):
            return (
                indexish(lineages, row["_lineage"]),
                indexish(isolates, row["mAb ID"]),
                row["mAb ID"])
        out.sort(key = sorter)
        for row in out:
            del row["_lineage"]
        with open(output[0], "w") as f_out:
            writer = csv.DictWriter(
                f_out,
                ["Animal ID", "Timepoint", "mAb ID", "Heavy Ig genes", "Light Ig genes",
                    "VH %nt mut", "VL %nt mut"],
                lineterminator="\n")
            writer.writeheader()
            writer.writerows(out)

def input_for_figS3(w):
    subject = {"A": "42056", "B": "6561", "C": "V031"}[w.panel]
    specimen = {"A": "RM42056WK72IGG", "B": "RM6561WK104IGG", "C": "RMV031WK56IGG"}[w.panel]
    targets = {
        "isolates": "metadata/isolates.csv",
        "junctions": f"from-paper/figS3{w.panel}_intermediate_junctions.txt",
        "summary_by_lineage": "analysis/isolates/summary_by_lineage.csv",
        "v": f"analysis/sonar/{subject}.IGH/{specimen}/germline/V.cysTruncated.fasta",
        "sonar": f"analysis/sonar/{subject}.IGH/{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa"}
    return targets

def params_for_figS3(w):
    lineage = {"A": "42056-a", "B": "6561-a", "C": "V031-a"}[w.panel]
    d_pos = {"A": 21, "B": 27, "C": 21}[w.panel]
    gap_pos = {"A": 0, "B": 49, "C": 34}[w.panel]
    gap_len = {"A": 0, "B": 3, "C": 6}[w.panel]
    return {"lineage": lineage, "d_pos": d_pos, "gap_pos": gap_pos, "gap_len": gap_len}

rule structure_output_figS3:
    input: expand("output/structure_figS3{panel}.fa", panel = ["A", "B", "C"])

rule structure_output_figS3panel:
    """FASTA for VDJ segments, intermediate lineage members, and mature isolated antibody for one lineage"""
    output: "output/structure_figS3{panel}.fa"
    input: unpack(input_for_figS3)
    params: params_for_figS3
    shell:
        """
            scripts/fig_s3.py {params[0][lineage]} {input.junctions} {params[0][d_pos]} {output} \
            --gap-pos {params[0][gap_pos]} --gap-len {params[0][gap_len]}
        """

rule structure_output_tableS3_d_info:
    """CSV of the raw germline D information that we summarized in table S3"""
    output: "output/structure_tableS3_d_info.csv"
    input: "analysis/aggregate-d/all.csv"
    run:
        subjects = [
            "5695", "6070", "6561", "40591", "41328",
            "42056", "44715", "T646", "V031", "V033"]
        with open(input[0]) as f_in, open(output[0], "w") as f_out:
            reader = csv.DictReader(f_in)
            rows = list(reader)
            rows.sort(key=lambda row: (subjects.index(row["subject"]), list(row.values())))
            writer = csv.DictWriter(f_out, reader.fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

### Info from structure paper itself

rule all_from_structure_paper:
    input: expand("from-structure-paper/{name}.csv", name=["fig1b", "fig2a", "tableS2"])

rule from_structure_paper:
    """Download CSV version of item from the structure paper from a copy in Google Sheets"""
    output: "from-structure-paper/{name}.csv"
    input: "metadata/paper_google_sheets.yml"
    shell: "scripts/download_google_sheet.py {wildcards.name} {input} {output}"
