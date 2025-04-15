### Overall workflow setup and metadata

import csv

# "germline-genbank" to use readymade IgDiscover outputs from GenBank, or
# "germline" to create them from scratch first
GERMLINE="germline-genbank"

# The Ramesh germline D sequences that we observed in the macaques that were
# not present in KIMDB.
RAMESH_D = ["IGHD1-1*01", "IGHD1-39*01", "IGHD3-26*02", "IGHD4-27*02"]

wildcard_constraints:
    # project metadata
    sample="[-A-Za-z0-9]+",
    specimen="[A-Za-z0-9]+",
    subject="[A-Za-z0-9]+",
    antibody_isolate=r"[-_A-Za-z0-9\.]+",
    antibody_lineage=r"[-_A-Za-z0-9\.]+",
    # antibody concepts
    chain="(heavy|light)",
    celltype="igm|igg",
    chain_type="(alpha|delta|gamma|mu|epsilon|kappa|lambda)",
    locus="(IGH|IGK|IGL)",
    segment="(V|D|J)",
    antibody_type="(IgA|IgD|IgG|IgM|IgE)",
    # other flow control
    name="[^/]+", # any match limited to one directory
    accession="[^/]+", # ditto, for accession numbers specifically
    word="[A-Za-z]+", # just letters

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

# from `sha256sum filenames`
def load_checksums():
    with open("metadata/sha256.txt") as f_in:
        pairs = dict(re.fullmatch(r"([^ ]+) +(.*)", x.strip()).groups()[::-1] for x in f_in)
    return pairs

CHECKSUMS = load_checksums()

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
    shell: "scripts/tabulate_seqs.py -p '{params.pattern}' -x subject=5695 {input} -o {output}"

rule genbank_igdiscover:
    """The new IgDiscover sequences for all subjects"""
    output: "analysis/genbank/igdiscover.csv"
    input: "genbank-placeholders/igdiscover.txt.gz"
    params:
        pattern=r"Rhesus Macaque (?P<subject>.*) antibody germline sequence for "
            r"locus (?P<locus>IG[HKL]) (?P<segment>[VJ]) segment"
    shell: "scripts/tabulate_seqs.py -p '{params.pattern}' -f gb {input} -o {output}"

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
            r"post-infection, (?P<locus>IG[HKL]) sequence"
    shell: "scripts/tabulate_seqs.py -p '{params.pattern}' -f gb {input} -o {output}"

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
