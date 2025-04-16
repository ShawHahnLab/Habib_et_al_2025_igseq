### SONAR (Lineage tracing with IgG reads)

from math import floor, log10

def make_sonar_rules():
    """Set up helper SONAR rules for each subject"""
    subject_locus_specs = {}
    subject_locus_specs_igg = {}
    def gather(obj, row):
        subject = row["igseq_Specimen_Subject"]
        locus = {"kappa": "IGK", "lambda": "IGL"}.get(row["igseq_Type"], "IGH")
        if subject not in obj:
            obj[subject] = {}
        if locus not in obj[subject]:
            obj[subject][locus] = []
        obj[subject][locus].append(row["igseq_Specimen"])
    for row in METADATA["biosamples"]:
        gather(subject_locus_specs, row)
        if "IgG" in row["igseq_Specimen_CellType"]:
            gather(subject_locus_specs_igg, row)
    for subject in sorted(subject_locus_specs):
        # annotation of all reads per biological specimen per locus, for both
        # IgG+ and IgM+ material
        targets = []
        for locus in ("IGH", "IGK", "IGL"):
            specimens = subject_locus_specs[subject].get(locus, [])
            targets += expand(
                "analysis/sonar/{subject}.{locus}/{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
                subject=subject, locus=locus, specimen=specimens)
        rule:
            name: f"sonar_module_1_{subject}"
            input: targets
        # Just IgG+, ID/DIV tables for comparing repertoire reads to germline V
        # and known antibodies of interest for whatever lineage(s) per subject
        targets = []
        for locus in ("IGH", "IGK", "IGL"):
            specimens = subject_locus_specs_igg[subject].get(locus, [])
            targets += expand(
                "analysis/sonar/{subject}.{locus}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab",
                subject=subject, locus=locus, specimen=specimens)
        rule:
            name: f"sonar_module_2_id_div_{subject}"
            input: targets
    # The automatic alignment tree outputs for each lineage, both heavy and light
    for row in METADATA["lineages"]:
        lineage = row["antibody_lineage"]
        targets = expand(
            "analysis/sonar/{subject}.{locus}/longitudinal.auto.{antibody_lineage}/"
            "output/longitudinal.auto.{antibody_lineage}_igphyml.tree",
            subject=row["subject"],
            locus=["IGH", row["light_locus"]],
            antibody_lineage=lineage)
        rule:
            name: f"sonar_module_3_igphyml_{lineage}"
            input: targets

make_sonar_rules()

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

def input_for_sonar_input(w):
    samples = []
    for row in METADATA["biosamples"]:
        locus = {"kappa": "IGK", "lambda": "IGL"}.get(row["igseq_Type"], "IGH")
        if row["igseq_Specimen"] == w.specimen and locus == w.locus:
            samples.append(row["sample_name"])
    return expand("analysis/merge/{sample}.fastq.gz", sample=samples)

rule all_sonar_module_1:
    input: make_target_sonar()

rule all_sonar_module_2_id_div:
    input: make_target_sonar("analysis/sonar/{subject}.{locus}/{specimen}/output/tables/{specimen}_goodVJ_unique_id-div.tab")

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
        V=f"analysis/{config['germline']}/{{subject}}.{{locus}}/V.fasta",
        D=f"analysis/{config['germline']}/{{subject}}.{{locus}}/D.fasta",
        J=f"analysis/{config['germline']}/{{subject}}.{{locus}}/J.fasta",
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
    input: f"analysis/{config['germline']}/{{subject}}.{{locus}}/V.fasta"
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
    # gets all the isolate seqs for one lineage and save as:
    #  * a text file with the IDs, so they can be given to sonar get_island
    #    like "--mab ID1 --mab ID2 . .." (as opposed to "--mab =a" for all of
    #    them from the ID/DIV table)
    #  * a FASTA file, so they can be given to sonar igphyml like "--natives
    #    path.fasta"
    output: "analysis/sonar/{subject}.{locus}/{name}/mab/mab.{antibody_lineage}.{word}"
    run:
        seq_col = "heavy_sequence" if wildcards.locus == "IGH" else "light_sequence"
        seen = {""} # same logic as for sonar_gather_mature
        with open(output[0], "wt") as f_out:
            for attrs in METADATA["isolates"]:
                seq = attrs[seq_col]
                seqid = attrs["antibody_isolate"]
                if attrs["subject"] == wildcards.subject \
                        and attrs["antibody_lineage"] == wildcards.antibody_lineage:
                    if wildcards.word == "txt":
                        if seq not in seen:
                            f_out.write(f"{seqid}\n")
                            seen.add(seq)
                    elif wildcards.word == "fasta":
                        f_out.write(f">{seqid}\n{seq}\n")
                    else:
                        raise ValueError

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

def sonar_module_3_collect_inputs(w):
    """List the inputs needed for SONAR module 3's collection step.

    This gives a dictionary of timepoint -> islandSeqs.fa pairs.
    """
    # This got quite convoluted in the main igseqhelper repo, with edge cases
    # where we have multiple distinct physical specimens isolated at the same
    # timepoint, but here I think it can be much simpler, since we have exactly
    # one IgG specimen per subject per timepoint.
    #
    # timepoint integers will be formatted like:
    # -12   "wkN12" # only V033 has negatives though
    #  -7   "wkN07"
    #   0   "wk000"
    #   4   "wk004"
    #   8   "wk008"
    #  12   "wk012"
    timepoint_specs = {}
    for row in METADATA["biosamples"]:
        if "IgG" in row["igseq_Specimen_CellType"] and w.subject == row["igseq_Specimen_Subject"]:
            tp_num = int(row["igseq_Specimen_Timepoint"])
            if tp_num in timepoint_specs and timepoint_specs[tp_num] != row['igseq_Specimen']:
                raise ValueError(f"Multiple timepoints for {tp_num} for {row['igseq_Specimen']}")
            timepoint_specs[tp_num] = row["igseq_Specimen"]
    # how many digits do we need to pad them all evenly? include an extra for
    # "-" if any are negative
    padlen = max((tp<0) + floor(log10(max(1, abs(tp))) + 1) for tp in timepoint_specs)
    # week such and such, "N" for negatives
    labels = ["wk" + str(num).replace("-", "N").zfill(padlen) for num in timepoint_specs]
    # sort by time, keeping time/label/specimens matched up
    trios = sorted(zip(timepoint_specs.keys(), labels, timepoint_specs.values()))
    # give dictionary of label to path pairs
    targets = expand(
        "analysis/sonar/{subject}.{locus}/{specimen}/"
        "output/sequences/nucleotide/islandSeqs_{antibody_lineage}.fa",
        subject=w.subject,
        locus=w.locus,
        specimen=[trio[2] for trio in trios],
        antibody_lineage=w.antibody_lineage)
    targets = {label: target for label, target in zip((trio[1] for trio in trios), targets)}
    return targets

def sonar_module_3_collect_param_seqs(_, input):
    args = [" --labels {key} --seqs {val}".format(key=key, val=Path(val).resolve()) for key, val in input.items()]
    return " ".join(args)

rule sonar_module_3_collect:
    output:
        collected="analysis/sonar/{subject}.{locus}/longitudinal.{word}.{antibody_lineage}/output/sequences/nucleotide/longitudinal.{word}.{antibody_lineage}-collected.fa"
    input: unpack(sonar_module_3_collect_inputs)
    singularity: "docker://jesse08/sonar"
    params:
        wd_sonar=lambda w: expand("analysis/sonar/{subject}.{locus}/longitudinal.{word}.{antibody_lineage}", **w),
        seqs=sonar_module_3_collect_param_seqs
    shell:
        """
            cd {params.wd_sonar}
            sonar merge_time {params.seqs}
        """

def input_for_sonar_module_3_igphyml(w):
    segs = ["V", "D", "J"] if w.locus == "IGH" else ["V", "J"]
    proj = "longitudinal.{word}.{antibody_lineage}"
    sonar_dir = "analysis/sonar/{subject}.{locus}/" + proj + "/"
    if w.word == "auto":
        targets = {seg: f"analysis/{config['germline']}/{w.subject}.{w.locus}/{seg}.fasta" for seg in segs}
        targets.update({
            "collected": sonar_dir + "output/sequences/nucleotide/" + proj + "-collected.fa",
            "natives": sonar_dir + "mab/mab.{antibody_lineage}.fasta"
        })
    else:
        targets = {
            "alignment": "analysis/sonar/{subject}.{locus}/alignment.{antibody_lineage}.{word}.fa"
        }
    return targets

def sonar_module_3_igphyml_param_v_id(wildcards):
    if dict(wildcards).get("word") != "auto":
        return ""
    lineages = {row["antibody_lineage"]: row for row in METADATA["lineages"]}
    key = "vl" if wildcards.locus in ["IGK", "IGL"] else "vh"
    v_call = lineages[wildcards.antibody_lineage][key]
    return v_call

rule sonar_module_3_igphyml:
    """SONAR 3: Run IgPhyML with auto or custom alignment and generate tree across specimens.

    For the automatic alignment case ("auto" as the keyword before the lineage
    name), this will build a tree with the appropriate germline V sequence for
    the lineage included and will root the tree on that V sequence.

    For the custom alignment case (any other keyword), this will given the
    first sequence ID in the alignment as the --root for sonar igphyml.
    """
    output:
        tree="analysis/sonar/{subject}.{locus}/longitudinal.{word}.{antibody_lineage}/output/longitudinal.{word}.{antibody_lineage}_igphyml.tree",
        inferred_nucl="analysis/sonar/{subject}.{locus}/longitudinal.{word}.{antibody_lineage}/output/sequences/nucleotide/longitudinal.{word}.{antibody_lineage}_inferredAncestors.fa",
        inferred_prot="analysis/sonar/{subject}.{locus}/longitudinal.{word}.{antibody_lineage}/output/sequences/amino_acid/longitudinal.{word}.{antibody_lineage}_inferredAncestors.fa",
        stats="analysis/sonar/{subject}.{locus}/longitudinal.{word}.{antibody_lineage}/output/logs/longitudinal.{word}.{antibody_lineage}_igphyml_stats.txt",
        # alignment produced here if word=auto, just copied from input (so the output is consistent) otherwise
        alignment="analysis/sonar/{subject}.{locus}/longitudinal.{word}.{antibody_lineage}/work/phylo/longitudinal.{word}.{antibody_lineage}_aligned.afa"
    input: unpack(input_for_sonar_module_3_igphyml)
    log: "analysis/sonar/{subject}.{locus}/longitudinal.{word}.{antibody_lineage}/log.txt"
    singularity: "docker://jesse08/sonar"
    params:
        wd_sonar=lambda w: expand("analysis/sonar/{subject}.{locus}/longitudinal.{word}.{antibody_lineage}", **w),
        seed=config.get("sonar_igphyml_seed", 123),
        aln_type=lambda w: "auto" if w.word == "auto" else "custom",
        args="-f",
        input_germline_v=lambda w, input: Path(input.V).resolve() if w.word == "auto" else "",
        input_natives=lambda w, input: Path(input.natives).resolve() if w.word == "auto" else "",
        input_alignment=lambda w, input: Path(input.alignment).resolve() if w.word != "auto" else "",
        output_alignment=lambda w, output: Path(output.alignment).resolve(),
        v_id=sonar_module_3_igphyml_param_v_id,
    shell:
        """
            # the perl inside the singularity container crashes with some
            # locale-related error if we leave LANG set, but it seems fine with
            # it unset
            unset LANG
            (
              cd {params.wd_sonar}
              date
              echo "Running SONAR Module 3 igphyml"
              echo
              echo "Alignment type: {params.aln_type}"
              echo "$(which sonar): $(sonar --version)"
              echo "Project directory: $PWD"
              echo "Random seed: {params.seed}"
              if [[ {params.aln_type} == "auto" ]]; then
                echo "Germline V ID for root: {params.v_id}"
                echo "Germline V library: {params.input_germline_v}"
                echo "Natives FASTA: {params.input_natives}"
                set -x
                sonar igphyml \
                    -v '{params.v_id}' \
                    --lib {params.input_germline_v} \
                    --natives {params.input_natives} \
                    --seed {params.seed} \
                    {params.args} 2>&1
              else
                root=$(head -n 1 {params.input_alignment} | cut -c 2- | cut -f 1 -d ' ')
                echo "Detected seq ID for tree root from first FASTA record: $root"
                set -x
                sonar igphyml --root "$root" -i {params.input_alignment} --seed {params.seed} {params.args} 2>&1
                cp {params.input_alignment} {params.output_alignment}
              fi

            ) | tee -a {log}
        """
