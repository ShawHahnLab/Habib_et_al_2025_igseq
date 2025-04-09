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
    input: f"analysis/{GERMLINE}/{{subject}}.{{locus}}/V.fasta"
    conda: "igseq.yml"
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
