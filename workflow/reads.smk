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
        samples="analysis/trim/samples.csv"
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

rule trim_samples:
    output: "analysis/trim/samples.csv"
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
