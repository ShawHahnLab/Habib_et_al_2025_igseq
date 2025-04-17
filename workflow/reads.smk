### Basic read processing

TARGET_DEMUX = expand("analysis/demux-link/{sample}.{rp}.fastq.gz", sample=[row["sample_name"] for row in METADATA["biosamples"]], rp=["R1", "R2"])
TARGET_MERGE = expand("analysis/merge/{sample}.fastq.gz", sample=[row["sample_name"] for row in METADATA["biosamples"]])

rule all_demux:
    input: TARGET_DEMUX

rule all_merge:
    input: TARGET_MERGE

if config["rundata"] == "sra":
    # if "sra", we'll shortcut the first couple of processing steps (custom
    # bcl2fastq and demultiplexing) and download demultiplexed files from the
    # SRA
    TARGET_SRA = expand("analysis/from-sra/{srr}_{num}.fastq.gz", srr=[row["Run"] for row in METADATA["biosamples"]], num=[1, 2])
    rule all_sra:
        input: TARGET_SRA
else:
    # If anything other than literally "sra", take it to be a path where
    # Illumina run directories are stored, and define rules to use those files
    # for the initial processing.
    rule getreads:
        output:
            outdir=directory("analysis/reads/{run}"),
            r1="analysis/reads/{run}/Undetermined_S0_L001_R1_001.fastq.gz",
            i1="analysis/reads/{run}/Undetermined_S0_L001_I1_001.fastq.gz",
            r2="analysis/reads/{run}/Undetermined_S0_L001_R2_001.fastq.gz"
        input: ancient(f"{config['rundata']}/{{run}}")
        log:
            conda="analysis/reads/{run}/conda_build.txt"
        conda: "igseq.yml"
        threads: 28
        shell:
            """
                conda list --explicit > {log.conda}
                igseq getreads -t {threads} --threads-load $(({threads}<4 ? {threads} : 4)) {input}
            """
    # (demux rules are separate for each run, since the output files vary from
    # one run to the next)
    by_run = defaultdict(list)
    for row in METADATA["biosamples"]:
        by_run[row["igseq_Run"]].append(row)
    for runid, rows in by_run.items():
        samp_names = [row["sample_name"] for row in rows]
        rule:
            name: f"demux_{runid}"
            output:
                outdir=directory(f"analysis/demux/{runid}"),
                fqgz=expand(
                    "analysis/demux/{run}/{samp}.{rp}.fastq.gz",
                    run=runid, samp=samp_names, rp=["R1", "R2", "I1"])
            input:
                reads=expand(
                    "analysis/reads/{run}/Undetermined_S0_L001_{rp}_001.fastq.gz",
                    run=runid, rp=["I1", "R1", "R2"]),
                samples=ancient("analysis/igseq_samples.csv")
            log:
                conda=f"analysis/demux/{runid}/conda_build.txt"
            conda: "igseq.yml"
            shell:
                """
                    conda list --explicit > {log.conda}
                    igseq demux --samples {input.samples} {input.reads}
                """

rule sra:
    output: expand("analysis/from-sra/{{srr}}_{rp}.fastq.gz", rp=[1, 2])
    shell: "fastq-dump --split-files --gzip {wildcards.srr} --outdir analysis/from-sra"

def input_for_demux_sample_link(w):
    for row in METADATA["biosamples"]:
        if row["sample_name"] == w.sample:
            break
    else:
        raise ValueError(w.sample)
    if config["rundata"] != "sra":
        runid = row["igseq_Run"]
        return {f"r{p}": f"analysis/demux/{runid}/{w.sample}.R{p}.fastq.gz" for p in (1, 2)}
    else:
        sra_Run = row["Run"]
        return {f"r{p}": f"analysis/from-sra/{accession}.R{p}.fastq.gz" for p in (1, 2)}

rule demux_link:
    # symlink files to analysis/demux/{sample} via either locally-run
    # demultiplexing (if using raw Illumina run data) or SRA-downloaded files
    output:
        r1="analysis/demux-link/{sample}.R1.fastq.gz",
        r2="analysis/demux-link/{sample}.R2.fastq.gz"
    input: unpack(input_for_demux_sample_link)
    shell:
        """
            ln -sr {input.r1} {output.r1}
            ln -sr {input.r2} {output.r2}
        """

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
        samples="analysis/igseq_samples.csv"
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

rule igseq_samples:
    output: "analysis/igseq_samples.csv"
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
