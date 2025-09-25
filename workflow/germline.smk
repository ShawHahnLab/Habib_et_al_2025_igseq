### IgDiscover (individualized germline references with IgM reads)

def make_germline_rules():
    """Set up helper IgDiscover and MINING-D rules for each subject"""
    subject_loci = defaultdict(set)
    for row in METADATA["biosamples"]:
        if "IgM" in row["igseq_Specimen_CellType"]:
            subject = row["igseq_Specimen_Subject"]
            locus = {"kappa": "IGK", "lambda": "IGL"}.get(row["igseq_Type"], "IGH")
            subject_loci[subject].add(locus)
    for subject, loci in subject_loci.items():
        targets_germline = expand("analysis/germline/{subject}.{locus}/{segment}.fasta",
            subject=subject, locus=sorted(loci), segment=("V", "D", "J"))
        rule:
            name: f"germline_{subject}"
            input: targets_germline
        targets_miningd = expand("analysis/mining-d/{subject}.output.{pval}.txt",
            subject=subject, pval=("default", "sensitive"))
        rule:
            name: f"miningd_{subject}"
            input: targets_miningd
        rule:
            name: f"germline_etc_{subject}"
            input: targets_germline + targets_miningd
    subject = "5695"
    targets_germline = expand("analysis/germline/{subject}.{locus}/{segment}.fasta",
        subject=subject, locus=("IGH", "IGL"), segment=("V", "D", "J"))
    rule:
        name: f"germline_{subject}"
        input: targets_germline

make_germline_rules()

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
        ["analysis/igdiscover/{ref}/{locus}/{subject}/stats/stats.json",
        "analysis/igdiscover/{ref}/{locus}/{subject}/custom_j_discovery/J.fasta"],
        zip, **attrs)

def make_target_miningd():
    subjects = set()
    for row in METADATA["biosamples"]:
        if row["igseq_Specimen_CellType"] == "IgM+":
            subjects.add(row["igseq_Specimen_Subject"])
    subjects = list(subjects)
    subjects.sort()
    return expand(
        "analysis/mining-d/{subject}.output.{pval}.txt",
        subject=subjects, pval=("default", "sensitive"))

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

TARGET_IGDISCOVER = make_target_igdiscover()
TARGET_MININGD = make_target_miningd()
TARGET_GERMLINE = make_target_germline()
TARGET_GERMLINE_GENBANK = make_target_germline(
    "analysis/germline-genbank/{subject}.{locus}/{{segment}}.fasta",
    (("IGH", "5695"), ("IGL", "5695")))

rule all_igdiscover:
    input: TARGET_IGDISCOVER

rule all_miningd:
    input: TARGET_MININGD

rule all_germline:
    input: TARGET_GERMLINE

rule all_germline_genbank:
    input: TARGET_GERMLINE_GENBANK

rule igdiscover_db_sonarramesh:
    """Prep IgDiscover starting database (SONAR/Ramesh for light chain)"""
    output: expand("analysis/igdiscover/sonarramesh/{{locus}}/{segment}.fasta", segment=["V", "D", "J"])
    params:
        outdir="analysis/igdiscover/sonarramesh/{locus}"
    conda: "igseq.yml"
    # (The grep process will keep only those starting with "IG" which will
    # exclude those prefixed with "ORF".  This way we don't need to manually
    # review and remove a spurious IGLV7-ADJ*01 from 40591's IGL
    # results at the end that appears to be expressed at low levels but
    # nonfunctional.)
    shell:
        """
            igseq show sonarramesh/{wildcards.locus}/{wildcards.locus}V | grep -A1 --no-group-sep '^>IG' > {output[0]}
            igseq show sonarramesh/IGH/IGHD > {output[1]}
            igseq show sonarramesh/{wildcards.locus}/{wildcards.locus}J > {output[2]}
        """

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
                row["Sample Name"]))
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
        db_j="analysis/igdiscover/{ref}/{locus}/{subject}/final/database/J.fasta",
        exp_d="analysis/igdiscover/{ref}/{locus}/{subject}/final/expressed_D.tab"
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

### MINING-D (D gene inference with IgM reads' CDR3s)

def input_for_miningd_get_cdrh3s(w):
    """List of SONAR AIRR tables for all IgM samples for one subject"""
    specs = []
    for attrs in METADATA["biosamples"]:
        if attrs["igseq_Specimen_Subject"] == w.subject and "IgM" in attrs["igseq_Specimen_CellType"]:
            specs.append(attrs["igseq_Specimen"])
    if not specs:
        raise ValueError(f"No specimens for MINING-D for {w.subject}?")
    return expand(
        "analysis/sonar/{subject}.IGH/{specimen}/output/tables/{specimen}_rearrangements.tsv",
        subject=w.subject, specimen=specs)

rule miningd_get_cdrh3s:
    """Gather CDRH3 sequences from SONAR's clustered antibody sequences"""
    # This uses SONAR simply to provide CDR3s of cleaned-up, deduplicated, and
    # clustered antibody sequences, even though we otherwise aren't generally
    # analyzing IgM+ material with SONAR.
    output: "analysis/mining-d/{subject}.cdr3.fasta"
    input: input_for_miningd_get_cdrh3s
    shell: "scripts/cdr3s_from_airr.py {input} {output}"

def checksum_for_miningd_run(w):
    out = expand("analysis/mining-d/{subject}.output.{pval}.txt", **dict(w))
    checksum = CHECKSUMS[out[0]]
    return checksum

# Note that the order and naming of sequences changes from run to run, but the
# sequences themselves do appear to typically be the same for any given input.
# Though, "typically" is the operative word; to handle the edge cases where
# there is some slight variation in the output, I'm using file hashes and
# Snakemake's retries feature to try to ensure we get the exact expected
# output.  (I *think* this is a result of MINING-D's use of Python sets and
# dictionaries but I'm not certain.)
rule miningd_run:
    """FASTA of inferred D sequences from MINING-D"""
    output:
        fasta="analysis/mining-d/{subject}.output.{pval}.fasta",
        txt=ensure("analysis/mining-d/{subject}.output.{pval}.txt", sha256=checksum_for_miningd_run)
    input: "analysis/mining-d/{subject}.cdr3.fasta"
    log:
        main="analysis/mining-d/{subject}.{pval}.log.txt",
        conda="analysis/mining-d/{subject}.{pval}.conda_build.txt"
    conda: "mining-d.yml"
    retries: 10
    threads: 8
    params:
        # The authors thought 600 was a good choice for rhesus macaque and
        # human
        num_k_mers=600,
        # 4.5e-36 is the default setting.  Will tend to leave off a few
        # nucleotides at the edges, but will avoid false positives.
        # (If a wildcard value is given that's not "default" or "sensitive",
        # take that text itself as the P value threshold.)
        p_val_th=lambda w: {"default": "4.5e-36", "sensitive": "4.5e-16"}.get(w.pval, w.pval)
    shell:
        """
            conda list --explicit > {log.conda}
            (
            numseqs_in=$(grep -c "^>" {input} || echo 0)
            echo "start: $(date +'%Y-%m-%d %H:%M:%S %Z')"
            echo "threads: {threads}"
            echo "input: {input}"
            echo "input seqs: $numseqs_in"
            echo "num_k_mers: {params.num_k_mers}"
            echo "p_val_th: {params.p_val_th}"
            MINING_D.py -t {threads} \
                -n {params.num_k_mers} \
                -p {params.p_val_th} \
                -i {input} -o {output.fasta}
            numseqs_out=$(grep -c "^>" {output.fasta} || echo 0)
            echo "output seqs: $numseqs_out"
            scripts/seq_list.py {output.fasta} {output.txt}
            echo "end: $(date +'%Y-%m-%d %H:%M:%S %Z')"
            ) >> {log.main}
        """


### Aggregated germline D information (for summarizing evidence across animals)

rule aggr_d_ref:
    """Make a D reference FASTA with KIMDB + select Ramesh entries"""
    output: "analysis/aggregate-d/d_ref.fasta"
    conda: "igseq.yml"
    params:
        ramesh_d_pat="(" + "|".join([x.replace("*", r"\*") for x in RAMESH_D]) + ")"
    shell:
        """
            igseq show kimdb/IGH/IGHD > {output}
            igseq show sonarramesh/IGH/IGHD | \
                grep -E -A1 --no-group-sep '{params.ramesh_d_pat}' | \
                sed -r 's/(>.*)/\\1_Ramesh/' >> {output}
        """

rule aggr_d:
    """Aggregate germline D information for one animal"""
    output: "analysis/aggregate-d/{subject}.csv"
    input:
        ref_d="analysis/aggregate-d/d_ref.fasta",
        igdisc_d="analysis/igdiscover/kimdb/IGH/{subject}/final/expressed_D.tab",
        min_d="analysis/mining-d/{subject}.output.sensitive.fasta"
    shell: "scripts/aggregate_d.py {input.ref_d} {input.igdisc_d} {input.min_d} {output}"

rule aggr_d_combined:
    output: "analysis/aggregate-d/all.csv"
    input: expand("analysis/aggregate-d/{subject}.csv", subject=sorted({attrs["igseq_Specimen_Subject"] for attrs in METADATA["biosamples"] if attrs["igseq_Specimen_CellType"] == "IgM+"}))
    params:
        subjects=sorted({attrs["igseq_Specimen_Subject"] for attrs in METADATA["biosamples"]})
    run:
        writer = None
        with open(output[0], "w") as f_out:
            for subject, path in zip(params.subjects, input):
                with open(path) as f_in:
                    reader = csv.DictReader(f_in)
                    if not writer:
                        writer = csv.DictWriter(
                            f_out,
                            fieldnames=["subject"] + list(reader.fieldnames),
                            lineterminator="\n")
                        writer.writeheader()
                    for row in reader:
                        row["subject"] = subject
                        writer.writerow(row)

### Germline (final versions of individualized germline references)

def input_for_germline(w):
    targets = {}
    ref="kimdb" if w.locus == "IGH" else "sonarramesh"
    for seg in ["V", "D", "J"]:
        targets[seg] = f"analysis/igdiscover/{ref}/{w.locus}/{w.subject}/final/database/{seg}.fasta"
    targets["filt_J"] = f"analysis/igdiscover/{ref}/{w.locus}/{w.subject}/custom_j_discovery/J.fasta"
    targets["ref_V"] = f"analysis/igdiscover/{ref}/{w.locus}/V.fasta"
    return targets

rule germline:
    """Prep germline files from IgDiscover's output"""
    output:
        V="analysis/germline/{subject}.{locus}/V.fasta",
        D="analysis/germline/{subject}.{locus}/D.fasta",
        J="analysis/germline/{subject}.{locus}/J.fasta"
    input: unpack(input_for_germline)
    # Special cases for 6561:
    #
    # For V, I spotted a "novel" germline sequence for IGL, IGLV4-ACF*02_S2122,
    # that's actually identical to to IGLV4-ACF*02 in the starting database.
    # It's not the only convoluted naming case like this, but I noticed this
    # one because lineage 6561-a uses it, so I reverted the sequence ID for
    # this one specifically.
    #
    # For J, 6561 is the only case where a J sequence in the starting DB is
    # subsequently filtered out by the custom_j_discovery logic.  Manual
    # inspection shows that 6561 evidently does really have a variant of IGHJ1
    # (KIMDB has IGHJ1-1*01 while Ramesh/IMGT has that same one as IGHJ1*01,
    # plus a synonymous variant IGHJ1*02 ending in TCCG versus TCAG) so I'm
    # bypassing the J filtering for 6561.
    shell:
        """
            if [[ {wildcards.subject} == 6561 ]]; then
                scripts/revert_igdiscover_names.py {input.ref_V} {input.V} {output.V}
                cp {input.J} {output.J}
            else
                cp {input.V} {output.V}
                cp {input.filt_J} {output.J}
            fi
            cp {input.D} {output.D}
        """

rule germline_link_5695:
    """Always use GenBank supplied IgDiscover files for 5695"""
    output:
        HV="analysis/germline/5695.IGH/V.fasta",
        HD="analysis/germline/5695.IGH/D.fasta",
        HJ="analysis/germline/5695.IGH/J.fasta",
        LV="analysis/germline/5695.IGL/V.fasta",
        LD="analysis/germline/5695.IGL/D.fasta",
        LJ="analysis/germline/5695.IGL/J.fasta"
    input:
        HV="analysis/germline-genbank/5695.IGH/V.fasta",
        HD="analysis/germline-genbank/5695.IGH/D.fasta",
        HJ="analysis/germline-genbank/5695.IGH/J.fasta",
        LV="analysis/germline-genbank/5695.IGL/V.fasta",
        LD="analysis/germline-genbank/5695.IGL/D.fasta",
        LJ="analysis/germline-genbank/5695.IGL/J.fasta"
    run:
        for path_in, path_out in zip(input, output):
            shell("ln -sr {path_in} {path_out}")

rule germline_genbank:
    """Prep germline files with the existing V and J from our GenBank entries"""
    output:
        V="analysis/germline-genbank/{subject}.{locus}/V.fasta",
        D="analysis/germline-genbank/{subject}.{locus}/D.fasta",
        J="analysis/germline-genbank/{subject}.{locus}/J.fasta"
    input:
        VJ="metadata/igdiscover.csv",
        # (In the light chain case the D file is just a placeholder to match
        # IgDiscover and provide a complete reference for igblast)
        D=lambda w: "analysis/igdiscover/" + ("kimdb" if w.locus == "IGH" else "sonarramesh") + "/{locus}/D.fasta",
        HJ="analysis/igdiscover/kimdb/IGH/J.fasta",
        LJ="analysis/igdiscover/sonarramesh/IGL/J.fasta"
    run:
        from Bio import SeqIO
        with open(input.VJ) as f_in, open(output.V, "w") as V_out, open(output.J, "w") as J_out:
            for row in csv.DictReader(f_in):
                if row["subject"] == wildcards.subject and row["locus"] == wildcards.locus:
                    handle = V_out if row["segment"] == "V" else J_out
                    seqid = row["sequence_id"]
                    seq = row["sequence"]
                    handle.write(f">{seqid}\n{seq}\n")
        seen = set()
        with open(input.D) as f_in, open(output.D, "w") as f_out:
            for rec in SeqIO.parse(f_in, "fasta"):
                if str(rec.seq) not in seen:
                    SeqIO.write(rec, f_out, "fasta-2line")
                    seen.add(str(rec.seq))
        # For 5695, use J from existing DB to get a complete reference.
        # (5695 is a special case since it was previously published so we
        # didn't re-deposit the reference in GenBank with these.)
        if wildcards.subject == "5695":
            ref = input.HJ if wildcards.locus == "IGH" else input.LJ
            shell("cp {ref} {output.J}")
