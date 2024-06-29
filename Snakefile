# See https://github.com/shawhahnlab/igseqhelper

import csv
from collections import defaultdict

# "germline-genbank" to use readymade IgDiscover outputs from GenBank, or
# "germline" to create them from scratch first
GERMLINE="germline-genbank"

wildcard_constraints:
    sample="[-A-Za-z0-9]+",
    specimen="[A-Za-z0-9]+",
    subject="[A-Za-z0-9]+",
    chain="(heavy|light)",
    celltype="igm|igg",
    chain_type="(alpha|delta|gamma|mu|epsilon|kappa|lambda)",
    locus="(IGH|IGK|IGL)",
    segment="(V|D|J)",
    antibody_type="(IgA|IgD|IgG|IgM|IgE)",
    antibody_isolate=r"[-_A-Za-z0-9\.]+",
    antibody_lineage=r"[-_A-Za-z0-9\.]+",
    accession="[^/]+",
    name="[^/]+",

### Setup

def load_metadata():
    metadata = {}
    for path in Path("metadata").glob("*.[tc]sv"):
        key = path.stem
        kwargs = {}
        if path.suffix == ".tsv":
            kwargs["delimiter"] = "\t"
        with open(path) as f_in:
            metadata[key] = list(csv.DictReader(f_in, **kwargs))
    for path in Path("metadata").glob("*.txt"):
        key = path.stem
        with open(path) as f_in:
            metadata[key] = [line.strip() for line in f_in]
    return metadata

METADATA = load_metadata()

def make_target_igdiscover():
    attrs = set()
    for row in METADATA["biosamples"]:
        if row["igseq_Specimen_CellType"] == "IgM+":
            ref = "kimdb" if row["igseq_Chain"] == "heavy" else "sonarramesh"
            locus = {"kappa": "IGK", "lambda": "IGL"}.get(row["igseq_Type"], "IGH")
            subject = row["igseq_Specimen_Subject"]
            attrs.add((ref, locus, subject))
    attrs = list(attrs)
    attrs.sort()
    attrs = {key: val for key, val in zip(("ref", "locus", "subject"), zip(*attrs))}
    return expand(
        "analysis/igdiscover/{ref}/{locus}/{subject}/stats/stats.json", zip, **attrs)

def make_target_germline(pattern="analysis/germline/{subject}.{locus}/{{segment}}.fasta", extra=None):
    attrs = set()
    for row in METADATA["biosamples"]:
        if row["igseq_Specimen_CellType"] == "IgM+":
            locus = {"kappa": "IGK", "lambda": "IGL"}.get(row["igseq_Type"], "IGH")
            attrs.add((locus, row["igseq_Specimen_Subject"]))
    attrs = list(attrs)
    attrs.sort()
    if extra:
        attrs.extend(extra)
    attrs = {key: val for key, val in zip(("locus", "subject"), zip(*attrs))}
    return expand(expand(pattern, zip, **attrs), segment=["V", "D", "J"])

def make_target_sonar(pattern="analysis/sonar/{subject}.{locus}/{specimen}/output/tables/{specimen}_rearrangements.tsv"):
    attrs = set()
    for row in METADATA["biosamples"]:
        if row["igseq_Specimen_CellType"] == "IgG+":
            locus = {"kappa": "IGK", "lambda": "IGL"}.get(row["igseq_Type"], "IGH")
            subject = row["igseq_Specimen_Subject"]
            specimen = row["igseq_Specimen"]
            attrs.add((subject, locus, specimen))
    attrs = list(attrs)
    attrs.sort()
    attrs = {key: val for key, val in zip(("subject", "locus", "specimen"), zip(*attrs))}
    return expand(pattern, zip, **attrs)

TARGET_TRIM = expand("analysis/trim/{sample}.fastq.gz", sample=[row["Sample"] for row in METADATA["samples"]])
TARGET_MERGE = expand("analysis/merge/{sample}.fastq.gz", sample=[row["Sample"] for row in METADATA["samples"]])
TARGET_IGDISCOVER = make_target_igdiscover()
TARGET_GERMLINE = make_target_germline()
TARGET_GERMLINE_GENBANK = make_target_germline(
    "analysis/germline-genbank/{subject}.{locus}/{{segment}}.fasta",
    (("IGH", "5695"), ("IGL", "5695")))
TARGET_SONAR_1 = make_target_sonar()
TARGET_SONAR_2_ID_DIV = make_target_sonar(
    "analysis/sonar/{subject}.{locus}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab")

rule all_sonar_1:
    input: TARGET_SONAR_1

rule all_sonar_2_id_div:
    input: TARGET_SONAR_2_ID_DIV

rule all_igdiscover:
    input: TARGET_IGDISCOVER

rule all_germline:
    input: TARGET_GERMLINE

rule all_germline_genbank:
    input: TARGET_GERMLINE_GENBANK

rule all_trim:
    input: TARGET_TRIM

rule all_merge:
    input: TARGET_MERGE

rule metadata_samples:
    output: "metadata/samples.csv"
    input: "metadata/biosamples.tsv"
    run:
        with open(input[0]) as f_in, open(output[0], "w") as f_out:
            writer = csv.DictWriter(
                f_out,
                ["Sample", "Run", "Specimen", "BarcodeFwd", "BarcodeRev", "Chain", "Type"],
                lineterminator="\n")
            writer.writeheader()
            for row in csv.DictReader(f_in, delimiter="\t"):
                row_out = {}
                row["igseq_Sample"] = row["sample_name"]
                for key, val in row.items():
                    key2 = key.removeprefix("igseq_")
                    if key2 in writer.fieldnames:
                        row_out[key2] = val
                writer.writerow(row_out)

rule download_genbank:
    output: "analysis/genbank/{accession}.{rettype}"
    shell: "scripts/download_ncbi.py nucleotide {wildcards.accession} {output} {wildcards.rettype}"

rule genbank_igdiscover_5695:
    """The old IgDiscover sequences for 5695 from the RHA1 paper

    We have a newer reference for 5695 now based on KIMDB for heavy chain, but
    here's what was used for the original RHA1 paper.
    """
    output: "analysis/genbank/igdiscover_5695.csv"
    input: expand("analysis/genbank/{accession}.fasta", accession=METADATA["genbank_igdiscover_5695"])
    params:
        pattern=r".*mulatta clone (?P<sequence_id>[^ ]+) immunoglobulin "
            r".*\((?P<locus>IG[HKL])(?P<segment>[VJ])\) mRNA"
    shell: "./scripts/tabulate_seqs.py -p '{params.pattern}' -x subject=5695 {input} -o {output}"

rule genbank_igdiscover:
    """The new IgDiscover sequences for all subjects"""
    output: "analysis/genbank/igdiscover.csv"
    input: "genbank-placeholders/igdiscover.txt.gz"
    params:
        pattern=r"Rhesus Macaque (?P<subject>.*) antibody germline sequence for "
            r"locus (?P<locus>IG[HKL]) (?P<segment>[VJ]) segment",
    shell: "./scripts/tabulate_seqs.py -p '{params.pattern}' -f gb {input} -o {output}"

rule genbank_isolates_5695:
    """Paired heavy and light chain sequences from 5695 from the RHA1 paper"""
    output: "analysis/genbank/isolates_5695.csv"
    input: expand("analysis/genbank/{accession}.fasta", accession=METADATA["genbank_isolates_5695"])
    # "We used an unbiased FACS strategy to isolate 20,000 individual memory B
    # cells from peripheral blood mononuclear cells (PBMCs) 65 weeks after SHIV
    # infection"
    params:
        pattern=r".*mulatta isolate (?P<antibody_isolate>[^ ]+) .* (?P<chain>[^ ]+) chain"
    shell:
        """
            ./scripts/tabulate_seqs.py -p '{params.pattern}' \
                -x subject=5695 -x antibody_lineage=5695-a -x timepoint=65 \
                {input} -o {output}
        """

rule genbank_isolates:
    """Paired heavy and light chain sequences across lineages"""
    output: "analysis/genbank/isolates.csv"
    input: expand("genbank-placeholders/isolates_{chain}.txt.gz", chain=["heavy", "light"])
    params:
        pattern=r"Rhesus macaque (?P<subject>.*) antibody lineage (?P<antibody_lineage>.*) "
            r"antibody (?P<antibody_isolate>.*) isolated at (?P<timepoint>[0-9]+) weeks "
            r"post-infection, (?P<locus>IG[HKL]) sequence",
    shell: "./scripts/tabulate_seqs.py -p '{params.pattern}' -f gb {input} -o {output}"

rule metadata_isolates:
    """Table of paired sequences for all lineages"""
    output: "metadata/isolates.csv"
    input:
        isolates="analysis/genbank/isolates.csv",
        isolates_5695="analysis/genbank/isolates_5695.csv"
    shell: "scripts/convert_gb_isolates.py {input} {output}"

rule metadata_igdiscover:
    """Table of new IgDiscover sequences for all subjects"""
    output: "metadata/igdiscover.csv"
    input: "analysis/genbank/igdiscover.csv"
    run:
        with open(output[0], "w") as f_out:
            writer = csv.DictWriter(
                f_out,
                ["subject", "locus", "segment", "sequence_id", "sequence"],
                lineterminator="\n")
            writer.writeheader()
            with open(input[0]) as f_in:
                for row in csv.DictReader(f_in):
                    row["sequence_id"] = row["sequence_id"].removeprefix(row["subject"] + "_")
                    writer.writerow(row)

### Basic read processing

rule trim:
    output:
        r1="analysis/trim/{sample}.R1.fastq.gz",
        r2="analysis/trim/{sample}.R2.fastq.gz",
        report1="analysis/trim/{sample}.cutadapt1.json",
        report2="analysis/trim/{sample}.cutadapt2.json",
        counts="analysis/trim/{sample}.trim.counts.csv"
    input:
        r1="analysis/demux/{sample}.R1.fastq.gz",
        r2="analysis/demux/{sample}.R2.fastq.gz",
        samples="metadata/samples.csv"
    log:
        conda="analysis/trim/{sample}.trim.conda_build.txt"
    conda: "igseq.yml"
    threads: 4
    shell:
        """
            conda list --explicit > {log.conda}
            igseq trim -t {threads} \
                --samples {input.samples} \
                --species rhesus \
                --outdir $(dirname {output.r1}) \
                {input.r1} {input.r2}
        """

rule merge:
    output:
        fqgz="analysis/merge/{sample}.fastq.gz",
        counts="analysis/merge/{sample}.merge.counts.csv",
        log="analysis/merge/{sample}.pear.log"
    input:
        r1="analysis/trim/{sample}.R1.fastq.gz",
        r2="analysis/trim/{sample}.R2.fastq.gz"
    log:
        conda="analysis/merge/{sample}.fastq.gz.conda_build.txt"
    conda: "igseq.yml"
    threads: 4
    shell:
        """
            conda list --explicit > {log.conda}
            igseq merge -t {threads} --outdir $(dirname {output.fqgz}) {input.r1} {input.r2}
        """

### IgDiscover (individualized germline references with IgM reads)

rule igdiscover_db_sonarramesh:
    """Prep IgDiscover starting database (SONAR/Ramesh for light chain)"""
    output: expand("analysis/igdiscover/sonarramesh/{{locus}}/{segment}.fasta", segment=["V", "D", "J"])
    params:
        outdir="analysis/igdiscover/sonarramesh/{locus}"
    conda: "igseq.yml"
    shell: "igseq vdj-gather sonarramesh/IGH/IGHD sonarramesh/{wildcards.locus} -o {params.outdir}"

rule igdiscover_db_kimdb:
    """Prep IgDiscover starting database (KIMDB 1.1 for heavy chain)"""
    output: expand("analysis/igdiscover/kimdb/IGH/{segment}.fasta", segment=["V", "D", "J"])
    conda: "igseq.yml"
    params: outdir="analysis/igdiscover/kimdb/IGH"
    shell: "igseq vdj-gather kimdb -o {params.outdir}"

def input_for_igdiscover_input(w):
    chain_type = {"IGH": "mu", "IGK": "kappa", "IGL": "lambda"}[w.locus]
    samples = []
    # Jumping through some hoops to ensure the ordering of samples matches what
    # I'd originally used in our analysis so that the input to IgDiscover is
    # byte-for-byte identical here as it was there.
    for row in METADATA["biosamples"]:
        if row["igseq_Specimen_CellType"] == "IgM+" and \
                row["igseq_Type"] == chain_type and \
                row["igseq_Specimen_Subject"] == w.subject:
            samples.append((
                METADATA["igm_runs"].index(row["igseq_Run"]),
                row["sample_name"]))
    samples = [s[1] for s in sorted(samples)]
    return expand("analysis/merge/{sample}.fastq.gz", sample=samples)

rule igdiscover_input:
    output: "analysis/igdiscover-input/{subject}.{locus}.fastq.gz"
    input: input_for_igdiscover_input
    run:
        if len(input) == 1:
            shell("ln -sr {input} {output}")
        elif input:
            shell("cat {input} > {output}")
        else:
            raise ValueError

rule igdiscover_init:
    output:
        yaml="analysis/igdiscover/{ref}/{locus}/{subject}/igdiscover.yaml",
        reads="analysis/igdiscover/{ref}/{locus}/{subject}/reads.fastq.gz",
    input:
        # IgDiscover always wants a "D" even for light
        db_v="analysis/igdiscover/{ref}/{locus}/V.fasta",
        db_d="analysis/igdiscover/{ref}/{locus}/D.fasta",
        db_j="analysis/igdiscover/{ref}/{locus}/J.fasta",
        reads="analysis/igdiscover-input/{subject}.{locus}.fastq.gz"
    conda: "igdiscover.yml"
    params:
        stranded="true",
        iterations=5
    shell:
        """
            rmdir $(dirname {output.yaml})
            igdiscover init \
                --db $(dirname {input.db_v}) \
                --single-reads {input.reads} \
                $(dirname {output.yaml})
            sed -i 's/^iterations: 1$/iterations: {params.iterations}/' {output.yaml}
            sed -i 's/^stranded: false$/stranded: {params.stranded}/' {output.yaml}
        """

rule igdiscover:
    output:
        stats="analysis/igdiscover/{ref}/{locus}/{subject}/stats/stats.json",
        db_v="analysis/igdiscover/{ref}/{locus}/{subject}/final/database/V.fasta",
        db_d="analysis/igdiscover/{ref}/{locus}/{subject}/final/database/D.fasta",
        db_j="analysis/igdiscover/{ref}/{locus}/{subject}/final/database/J.fasta"
    input:
        yaml="analysis/igdiscover/{ref}/{locus}/{subject}/igdiscover.yaml",
        r1="analysis/igdiscover/{ref}/{locus}/{subject}/reads.fastq.gz",
    log:
        conda="analysis/igdiscover/{ref}/{locus}/{subject}/igdiscover_run.conda_build.txt"
    conda: "igdiscover.yml"
    threads: 20
    shell:
        """
            conda list --explicit > {log.conda}
            cd $(dirname {output.stats})/.. && igdiscover run --cores {threads}
        """

rule igdiscover_custom_j_discovery:
    """Use more stringent approach for germline J inference"""
    output:
        tab="analysis/igdiscover/{ref}/{locus}/{subject}/custom_j_discovery/J.tab",
        fasta="analysis/igdiscover/{ref}/{locus}/{subject}/custom_j_discovery/J.fasta"
    input:
        # (This command actually uses iteration-01's J.fasta and
        # filtered.tsv.gz but I don't have those listed as part of the
        # input/output paths in all this.  So I'll just request the output of
        # my igdiscover rule.)
        db_j="analysis/igdiscover/{ref}/{locus}/{subject}/final/database/J.fasta",
    params:
        db_j="analysis/igdiscover/{ref}/{locus}/{subject}/iteration-01/database/J.fasta",
        tab="analysis/igdiscover/{ref}/{locus}/{subject}/iteration-01/filtered.tsv.gz",
        jcov=100,
        ratio=0.3
    log:
        conda="analysis/igdiscover/{ref}/{locus}/{subject}/custom_j_discovery/conda_build.txt"
    conda: "igdiscover.yml"
    shell:
        """
            arg_j_cov=""
            if [[ "{params.jcov}" != "" ]]; then
                arg_j_cov="--j-coverage {params.jcov}"
            fi
            arg_allele_ratio=""
            if [[ "{params.ratio}" != "" ]]; then
                arg_allele_ratio="--allele-ratio {params.ratio}"
            fi
            conda list --explicit > {log.conda}
            igdiscover discoverjd --database {params.db_j} \
                --gene J $arg_j_cov $arg_allele_ratio \
                {params.tab} --fasta {output.fasta} > {output.tab}
        """

### Germline (final versions of individualized germline references)

def input_for_germline(w):
    targets = {}
    ref="kimdb" if w.locus == "IGH" else "sonarramesh"
    for seg in ["V", "D", "J"]:
        pattern = "analysis/igdiscover/{ref}/{locus}/{subject}/" + \
            ("custom_j_discovery/J.fasta" if seg == "J" else "final/database/{seg}.fasta")
        targets[seg] = expand(pattern, ref=ref, locus=w.locus, subject=w.subject, seg=seg)
    return targets

rule germline:
    """Prep germline files from IgDiscover's output"""
    output:
        V="analysis/germline/{subject}.{locus}/V.fasta",
        D="analysis/germline/{subject}.{locus}/D.fasta",
        J="analysis/germline/{subject}.{locus}/J.fasta",
    input: unpack(input_for_germline)
    shell:
        """
            cp {input.V} {output.V}
            cp {input.D} {output.D}
            cp {input.J} {output.J}
        """

rule germline_genbank:
    """Prep germline files with the existing V and J from our GenBank entries"""
    output:
        V="analysis/germline-genbank/{subject}.{locus}/V.fasta",
        D="analysis/germline-genbank/{subject}.{locus}/D.fasta",
        J="analysis/germline-genbank/{subject}.{locus}/J.fasta",
    input:
        VJ="metadata/igdiscover.csv",
        D="analysis/igdiscover/kimdb/IGH/D.fasta",
        HJ="analysis/igdiscover/kimdb/IGH/J.fasta",
        LJ="analysis/igdiscover/sonarramesh/IGL/J.fasta"
    run:
        with open(input.VJ) as f_in, open(output.V, "w") as V_out, open(output.J, "w") as J_out:
            for row in csv.DictReader(f_in):
                if row["subject"] == wildcards.subject and row["locus"] == wildcards.locus:
                    handle = V_out if row["segment"] == "V" else J_out
                    seqid = row["sequence_id"]
                    seq = row["sequence"]
                    handle.write(f">{seqid}\n{seq}\n")
        shell("cp {input.D} {output.D}")
        # For 5695, use J from existing DB to get a complete reference.
        # (5695 is a special case since it was previously published so we
        # didn't re-deposit the reference in GenBank with these.)
        if wildcards.subject == "5695":
            ref = input.HJ if wildcards.locus == "IGH" else input.LJ
            shell("cp {ref} {output.J}")

### Isolates

def input_for_igblast_isolates_combined(w):
    attrs = set()
    for row in METADATA["isolates"]:
        attrs.add((row["subject"], "IGH"))
        if row["light_sequence"]:
            attrs.add((row["subject"], row["light_locus"]))
    attrs = list(attrs)
    attrs.sort()
    attrs = {key: val for key, val in zip(("subject", "locus"), zip(*attrs))}
    return expand("analysis/isolates/{subject}.{locus}/igblast.tsv", zip, **attrs)

rule igblast_isolates_combined:
    """Combined TSV of IgBLAST AIRR across heavy+light seqs of all isolates"""
    output: "analysis/isolates/igblast.tsv"
    input: input_for_igblast_isolates_combined
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
        ref=expand("analysis/{germline}/{{subject}}.{{locus}}/{segment}.fasta", germline=GERMLINE, segment=["V", "D", "J"])
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

### SONAR (Lineage tracing with IgG reads)

def input_for_sonar_input(w):
    samples = []
    for row in METADATA["biosamples"]:
        locus = {"kappa": "IGK", "lambda": "IGL"}.get(row["igseq_Type"], "IGH")
        if row["igseq_Specimen"] == w.specimen and locus == w.locus:
            samples.append(row["sample_name"])
    return expand("analysis/merge/{sample}.fastq.gz", sample=samples)

rule sonar_input:
    output: "analysis/sonar-input/{specimen}.{locus}.fastq.gz"
    input: input_for_sonar_input
    run:
        if len(input) == 1:
            shell("ln -sr {input} {output}")
        elif input:
            shell("cat {input} > {output}")
        else:
            raise ValueError

WD_SONAR = Path("analysis/sonar/{subject}.{locus}/{specimen}")

JMOTIF = {
    "IGH": "TGGGG",
    "IGK": "TT[C|T][G|A]G",
    "IGL": "TT[C|T][G|A]G"}

rule sonar_module_1:
    """SONAR module 1 (antibody sequence annotation and clustering)"""
    output:
        fasta=WD_SONAR/"output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        rearr=WD_SONAR/"output/tables/{specimen}_rearrangements.tsv"
    input:
        V=f"analysis/{GERMLINE}/{{subject}}.{{locus}}/V.fasta",
        D=f"analysis/{GERMLINE}/{{subject}}.{{locus}}/D.fasta",
        J=f"analysis/{GERMLINE}/{{subject}}.{{locus}}/J.fasta",
        reads="analysis/sonar-input/{specimen}.{locus}.fastq.gz"
    log: (WD_SONAR/"log.txt").resolve()
    singularity: "docker://jesse08/sonar"
    threads: 4
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        cluster_id_fract=.99,
        cluster_min2=2,
        libv_arg=lambda w, input: "--lib " + str(Path(input.V).resolve()),
        libd_arg=lambda w, input: "D" in input and "--dlib " + str(Path(input.D).resolve()) or "--noD",
        libj_arg=lambda w, input: "--jlib " + str(Path(input.J).resolve()),
        jmotif=lambda w: JMOTIF[w.locus]
    shell:
        """
            zcat {input.reads} > {params.wd_sonar}/reads.fastq
            cd {params.wd_sonar}
            date | tee -a {log}
            echo "$(which sonar): $(sonar --version)" | tee -a {log}
            echo "Running sonar module 1" | tee -a {log}
            echo "Project directory: $PWD" | tee -a {log}
            echo "germline V: {input.V}" | tee -a {log}
            echo "germline D: {input.D}" | tee -a {log}
            echo "germline J: {input.J}" | tee -a {log}
            echo "J motif: {params.jmotif}" | tee -a {log}
            echo "Cluster ID fract: {params.cluster_id_fract}" | tee -a {log}
            echo "Cluster min2: {params.cluster_min2}" | tee -a {log}
            echo | tee -a {log}
            sonar blast_V {params.libv_arg} --derep --threads {threads} 2>&1 | tee -a {log}
            sonar blast_J {params.libd_arg} {params.libj_arg} --noC --threads {threads} 2>&1 | tee -a {log}
            sonar finalize --jmotif '{params.jmotif}' --threads {threads} 2>&1 | tee -a {log}
            sonar cluster_sequences --id {params.cluster_id_fract} --min2 {params.cluster_min2} 2>&1 | tee -a {log}
        """

rule sonar_gather_mature:
    """Get heavy or light chain mature antibody sequences as references for SONAR."""
    output: WD_SONAR/"mab/mab.fasta"
    run:
        seq_col = "heavy_sequence" if wildcards.locus == "IGH" else "light_sequence"
        seen = {""} # skip duplicates below, but always skip empty entries too
        with open(output[0], "wt") as f_out:
            for attrs in METADATA["isolates"]:
                if seq_col == "light_sequence" and attrs["light_locus"] != wildcards.locus:
                    # Skip light chain sequences that are for the other locus
                    # than whatever was amplified here
                    continue
                seqid = attrs["antibody_isolate"]
                seq = attrs[seq_col]
                if attrs["subject"] == wildcards.subject and seq not in seen:
                    f_out.write(f">{seqid} {attrs['antibody_lineage']}\n")
                    f_out.write(attrs[seq_col]+"\n")
                    seen.add(seq)

rule sonar_v_trunc:
    """Truncate germline V sequences through the conserved C just before CDR3, for SONAR"""
    output: WD_SONAR/"germline/V.cysTruncated.fasta"
    input: f"analysis/{GERMLINE}/{{subject}}.{{locus}}/V.fasta",
    shell: "scripts/truncate_v.py {input} {output}"

rule sonar_module_2_id_div:
    """SONAR module 2 IDentity (to V) and DIVergence (to mature) calculations"""
    output: WD_SONAR/"output/tables/{specimen}_goodVJ_unique_id-div.tab"
    input:
        V=WD_SONAR/"germline/V.cysTruncated.fasta",
        fasta=WD_SONAR/"output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        mab=WD_SONAR/"mab/mab.fasta"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_v=lambda w, input: Path(input.V).resolve(),
        input_mab=lambda w, input: Path(input.mab).resolve()
    singularity: "docker://jesse08/sonar"
    threads: 4
    shell:
        """
            cd {params.wd_sonar}
            sonar id-div -g "{params.input_v}" -a "{params.input_mab}" -t {threads}
        """

rule sonar_list_members_for_lineage:
    # gets all the seq IDs as in the input mab antibody FASTA and save in a
    # text file, so they can be given to sonar get_island like "--mab ID1 --mab
    # ID2 ..." (as opposed to "--mab =a" for all of them from the ID/DIV table)
    output: WD_SONAR / "mab/mab.{antibody_lineage}.txt"
    run:
        seq_col = "heavy_sequence" if wildcards.locus == "IGH" else "light_sequence"
        seen = {""} # same logic as for sonar_gather_mature
        with open(output[0], "wt") as f_out:
            for attrs in METADATA["isolates"]:
                seq = attrs[seq_col]
                seqid = attrs["antibody_isolate"]
                if attrs["subject"] == wildcards.subject \
                        and seq not in seen and \
                        attrs["antibody_lineage"] == wildcards.antibody_lineage:
                    f_out.write(f"{seqid}\n")
                    seen.add(seq)

# NOTE this step is interactive over X11
##
# To get X working with Snakemake + singularity I had to append both of these
# arguments to the snakemake call (where our /home is mounted from /data/home):
#
#     --use-singularity --singularity-args "-B /home -B /data/home -H /home/$USER"
#
# I think these home directory options are needed so that the Xauthority stuff
# works but I'm not totally sure.   Using a regular Singularity image doesn't
# need this so I think it must be something about how Snakemake calls
# singularity.
rule sonar_module_2_id_div_island:
    output:
        seqids=WD_SONAR / "output/tables/islandSeqs_{antibody_lineage}.txt"
    input:
        iddiv=WD_SONAR / "output/tables/{specimen}_goodVJ_unique_id-div.tab",
        mab=WD_SONAR / "mab/mab.{antibody_lineage}.txt"
    singularity: "docker://jesse08/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_iddiv=lambda w, input: Path(input.iddiv).resolve(),
        outprefix=lambda w, output: Path(output.seqids).stem
    shell:
        """
            mabargs=$(sed 's/^/--mab /' {input.mab})
            cd {params.wd_sonar}
            sonar get_island {params.input_iddiv} $mabargs --output {params.outprefix} --plotmethod binned
        """

rule sonar_module_2_id_div_getfasta:
    """SONAR 2: Extract FASTA matching selected island's seq IDs."""
    output:
        fasta=WD_SONAR / "output/sequences/nucleotide/islandSeqs_{txt}.fa"
    input:
        fasta=WD_SONAR / "output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        seqids=WD_SONAR / "output/tables/islandSeqs_{txt}.txt"
    singularity: "docker://jesse08/sonar"
    params:
        wd_sonar=lambda w: expand(str(WD_SONAR), **w),
        input_fasta=lambda w, input: Path(input.fasta).resolve(),
        input_seqids=lambda w, input: Path(input.seqids).resolve(),
        output_fasta=lambda w, input, output: Path(output.fasta).resolve()
    shell:
        """
            cd {params.wd_sonar}
            sonar getFastaFromList \
                -l {params.input_seqids} \
                -f {params.input_fasta} \
                -o {params.output_fasta}
        """

### Output
#
# Some final summary outputs approximating what's shown in the paper itself.

def final_gene_name(txt):
    """Format the gene names like the paper has them"""
    return re.sub(r"IG([HKL])([VDJ])", r"\2\1", txt)

def indexish(items, item, default=-1):
    """Index of item in items or -1 if not present (just to help override sorting)"""
    try:
        return items.index(item)
    except ValueError:
        return default

rule output_fig1b:
    output: "output/fig1b.csv"
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

rule output_fig2a:
    output: "output/fig2a.csv"
    input: "analysis/isolates/summary.csv",
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

rule output_tableS2:
    output: "output/tableS2.csv"
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

rule output_figS3:
    input: expand("output/figS3{panel}.fa", panel = ["A", "B", "C"])

rule output_figS3A:
    output: "output/figS3{panel}.fa"
    input: unpack(input_for_figS3)
    params: params_for_figS3
    shell:
        """
            scripts/fig_s3.py {params[0][lineage]} {input.junctions} {params[0][d_pos]} {output} \
            --gap-pos {params[0][gap_pos]} --gap-len {params[0][gap_len]}
        """

### Info from paper itself

rule all_from_paper:
    input: expand("from-paper/{name}.csv", name=["fig1b", "fig2a", "tableS2"])

rule from_paper:
    """Download CSV version of item from the paper from a copy in Google Sheets"""
    output: "from-paper/{name}.csv"
    input: "metadata/paper_google_sheets.yml"
    shell: "scripts/download_google_sheet.py {wildcards.name} {input} {output}"
