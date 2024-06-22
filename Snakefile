# See https://github.com/shawhahnlab/igseqhelper

import csv

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
    return metadata

METADATA = load_metadata()

def make_target_igdiscover():
    attrs = set()
    for row in METADATA["biosamples"]:
        if row["igseq_Specimen_CellType"] == "IgM+":
            ref = "kimdb" if row["igseq_Chain"] == "heavy" else "sonarramesh"
            locus = {"mu": "IGH", "kappa": "IGK", "lambda": "IGL"}[row["igseq_Type"]]
            subject = row["igseq_Specimen_Subject"]
            attrs.add((ref, locus, subject))
    attrs = list(attrs)
    attrs.sort()
    attrs = {key: val for key, val in zip(("ref", "locus", "subject"), zip(*attrs))}
    return expand(
        "analysis/igdiscover/{ref}/{locus}/{subject}/stats/stats.json", zip, **attrs)

TARGET_IGDISCOVER = make_target_igdiscover()
TARGET_TRIM = expand("analysis/trim/{sample}.fastq.gz", sample=[row["Sample"] for row in METADATA["samples"]])
TARGET_MERGE = expand("analysis/merge/{sample}.fastq.gz", sample=[row["Sample"] for row in METADATA["samples"]])

rule all_igdiscover_run:
    input: TARGET_IGDISCOVER

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

### IgDiscover

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
    for row in METADATA["biosamples"]:
        if row["igseq_Specimen_CellType"] == "IgM+" and \
                row["igseq_Type"] == chain_type and \
                row["igseq_Specimen_Subject"] == w.subject:
            samples.append(row["sample_name"])
    return expand("analysis/merge/{sample}.fastq.gz", sample=samples)

rule igdiscover_input:
    output: "analysis/igdiscover-input/{subject}.{locus}.fastq.gz"
    input: input_for_igdiscover_input
    run:
        if len(input) == 1:
            Path(output[0]).symlink_to(Path(Path(input[0]).name)/inputs[0].name)
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

rule igdiscover_run:
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
