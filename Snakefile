# See https://github.com/shawhahnlab/igseqhelper

import csv

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
            locus = {"kappa": "IGK", "lambda": "IGL"}.get(row["igseq_Type"], "IGH")
            subject = row["igseq_Specimen_Subject"]
            attrs.add((ref, locus, subject))
    attrs = list(attrs)
    attrs.sort()
    attrs = {key: val for key, val in zip(("ref", "locus", "subject"), zip(*attrs))}
    return expand(
        "analysis/igdiscover/{ref}/{locus}/{subject}/stats/stats.json", zip, **attrs)

def make_target_sonar_1():
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
    return expand(
        "analysis/sonar/{subject}.{locus}/{specimen}/"
        "output/tables/{specimen}_rearrangements.tsv", zip, **attrs)

TARGET_TRIM = expand("analysis/trim/{sample}.fastq.gz", sample=[row["Sample"] for row in METADATA["samples"]])
TARGET_MERGE = expand("analysis/merge/{sample}.fastq.gz", sample=[row["Sample"] for row in METADATA["samples"]])
TARGET_IGDISCOVER = make_target_igdiscover()
TARGET_SONAR_1 = make_target_sonar_1()

rule all_sonar_1:
    input: TARGET_SONAR_1

rule all_igdiscover:
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

rule metadata_isolates:
    output: "metadata/isolates.csv"
    input: expand("genbank-placeholders/isolates_{chain}.txt.gz", chain=["heavy", "light"])
    shell: "scripts/convert_gb_isolates.py {input} {output}"

rule metadata_igdiscover:
    output: "metadata/igdiscover.csv"
    input: "genbank-placeholders/igdiscover.txt.gz"
    shell: "scripts/convert_gb_igdiscover.py {input} {output}"

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

rule germline:
    output:
        V="analysis/germline/{subject}.{locus}/V.fasta",
        D="analysis/germline/{subject}.{locus}/D.fasta",
        J="analysis/germline/{subject}.{locus}/J.fasta",
    input: unpack(lambda w: {seg: expand("analysis/igdiscover/{ref}/{locus}/{subject}/final/database/{seg}.fasta", ref="kimdb" if w.locus == "IGH" else "sonarramesh", locus=w.locus, subject=w.subject, seg=seg) for seg in ["V", "D", "J"]})
    shell:
        """
            cp {input.V} {output.V}
            cp {input.D} {output.D}
            cp {input.J} {output.J}
        """

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
            Path(output[0]).symlink_to(Path(Path(input[0]).name)/inputs[0].name)
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
    output:
        fasta=WD_SONAR/"output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        rearr=WD_SONAR/"output/tables/{specimen}_rearrangements.tsv"
    input:
        V="analysis/germline/{subject}.{locus}/V.fasta",
        D="analysis/germline/{subject}.{locus}/D.fasta",
        J="analysis/germline/{subject}.{locus}/J.fasta",
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
