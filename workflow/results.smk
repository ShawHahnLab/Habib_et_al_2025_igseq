### gathering final results in one place

TARGET_RESULTS_LINEAGES = expand(
    expand("results/{subject}/{antibody_lineage}/{antibody_lineage}_{{chain}}_igphyml.tree",
    zip,
    subject=[row["subject"] for row in METADATA["lineages"] if row["subject"] != "5695"],
    antibody_lineage=[row["antibody_lineage"] for row in METADATA["lineages"] if row["subject"] != "5695"]),
    chain=("heavy", "light"))

rule all_results:
    input: TARGET_RESULTS_LINEAGES

def make_results_rules():
    for row in METADATA["lineages"]:
        if row["subject"] == "5695":
            continue
        targets = expand("results/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_{thing}",
            subject=row["subject"],
            antibody_lineage=row["antibody_lineage"],
            chain=("heavy", "light"),
            thing=("aligned.afa", "igphyml.tree", "collected.fa", "inferredAncestors.fa"))
        rule:
            name: f"results_lineage_{row['antibody_lineage']}"
            input: targets

make_results_rules()

# Subject summary rules

# Antibody Lineage summary rules

def set_locus(w):
    """Set the locus wildcard based on chain and lineage and return as dict"""
    w = dict(w)
    for attrs in METADATA["lineages"]:
        if w["antibody_lineage"] == attrs["antibody_lineage"]:
            break
    else:
        raise ValueError(f"Can't find lineage {antibody_lineage}")
    w["locus"] = attrs["light_locus"] if w["chain"] != "heavy" else "IGH"
    return w

rule results_lineage_tree:
    output: "results/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_igphyml.tree"
    input: lambda w: expand("analysis/sonar/{subject}.{locus}/longitudinal.auto.{antibody_lineage}/output/longitudinal.auto.{antibody_lineage}_igphyml.tree", **set_locus(w))
    shell: "cp {input} {output}"

# Using `igseq convert` for FASTA files to ensure we get unwrapped versions
rule results_lineage_aln_auto:
    output: "results/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_aligned.afa"
    input: lambda w: expand("analysis/sonar/{subject}.{locus}/longitudinal.auto.{antibody_lineage}/work/phylo/longitudinal.auto.{antibody_lineage}_aligned.afa", **set_locus(w))
    conda: "igseq.yml"
    shell: "igseq convert {input} {output}"
    
rule results_lineage_collected:
    output: "results/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_collected.fa"
    input: lambda w: expand("analysis/sonar/{subject}.{locus}/longitudinal.auto.{antibody_lineage}/output/sequences/nucleotide/longitudinal.auto.{antibody_lineage}-collected.fa", **set_locus(w))
    conda: "igseq.yml"
    shell: "igseq convert {input} {output}"

rule results_lineage_ancestors_auto:
    output: "results/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_inferredAncestors.fa"
    input: lambda w: expand("analysis/sonar/{subject}.{locus}/longitudinal.auto.{antibody_lineage}/output/sequences/nucleotide/longitudinal.auto.{antibody_lineage}_inferredAncestors.fa", **set_locus(w))
    conda: "igseq.yml"
    shell: "igseq convert {input} {output}"
