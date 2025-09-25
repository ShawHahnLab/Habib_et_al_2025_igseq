# Habib et al 2025: V2 Apex Antibody Lineage Tracing

Analysis code and related metadata for the antibody lineage tracing aspect of
our 2025 HIV/V2 apex antibody coevolution paper.  The workflow is built with
[Snakemake] and [Python] to handle the steps of the analysis, [conda] for
software dependency management, and [Singularity] for a [Docker] container for
one aspect (SONAR).  This code makes use of [igseq] for some basic data
processing steps, which handles [bcl2fastq] and custom Python code for
specialized demultiplexing, [cutadapt] for adapter trimming, and [PEAR] for
read merging.  (The material in this repository is largely a streamlined
version of our [igseqhelper] workflow used for ongoing lineage tracing
analyses.)

Associated information from the publication:

 * Antibody repertoire NGS SRA project: PRJNA1121265 (BioSamples TSV provided
   here as metadata/biosamples.tsv)
 * Mature antibody paired heavy and light chain GenBank entries: MT610888 -
   MT610895, PP909817 - PP910010 (provided here as metadata/isolates.csv)
 * Per-subject antibody germline sequence reference GenBank entries: PP648252 -
   PP649987 (provided here as metadata/igdiscover.csv)
 * Antibody lineage UCA inferred heavy and light chain GenBank entries:
   PV467379 - PV467402

File layout here:

 * `metadata/`: lists and tables defining samples and sample metadata and
   related information
 * `results/`: one directory per subject containing germline-related analysis
   output files, with subdirectories per lineage containing lineage-specific
   outputs
 * `workflow/`: Snakemake workflow files (.smk) and associated conda
   environment definitions (.yml)
 * `scripts/`: Python scripts used in the workflow steps
 * `environment.yml`: a top-level conda environment definition for the
   repository
 * `analysis/`: the Snakemake rules for the workflow will write analysis files
   to subdirectories here
 * `structure-paper/`: text equivalents of a few items from our closely-related
   [V2 apex antibody structure paper], and recreated outputs with code here

Software usage, assuming Linux and an existing [mamba] install for conda package
management:

    $ git clone https://github.com/ShawHahnLab/Habib_et_al_2025_igseq.git
    $ cd Habib_et_al_2025_igseq
    $ mamba env update --file environment.yml # will create v2-apex-abs-paper environment
    $ mamba activate v2-apex-abs-paper
    $ snakemake --use-conda --use-singularity {output file paths or snakemake rules}

Overall workflow:

 1. Initial read processing (`workflow/reads.smk`)
 2. IgDiscover and MINING-D for germline sequence analysis; to use predefined
    germline files in downstream steps use `--config germline=germline-genbank`
    when calling snakemake
 3. SONAR for lineage tracing (`workflow/sonar.smk`)
    1. SONAR "module 1" (annotate) for read annotation and clustering, with one
       SONAR project directory per subject, locus, and timepoint
    2. SONAR "module 2" (lineage) for comparison of clustered and filtered
       repertoire sequences to known antibodies ("native" mAbs from the paired
       isolate sequences provided) for lineages of interest, creating
       per-lineage files in each SONAR project directory from module 1
    3. SONAR "module 3" (phylogeny) for longitudinal phylogenetic analysis
       combining lineage member sequences discovered in module 2 with the
       known "native" antibodies; creates a separate project directory per
       subject, locus, and lineage

Some key processing rules, named by subject or lineage:

    germline_{subject}
    sonar_module_1_{subject}
    sonar_module_2_id_div_{subject}
    sonar_module_3_igphyml_{lineage}

For example:

    snakemake -c 20 --use-conda --use-singularity sonar_module_2_id_div_6561

Or to run all steps necessary to produce heavy and light chain IgPhyML trees
and associated files for lineage 6561-a, with automatic alignment and default
(germline-V-based) tree rooting:

    snakemake -c 20 --use-conda --use-singularity sonar_module_3_igphyml_6561-a

Note that one step, selecting "islands" of lineage members from SONAR's ID/DIV
plots, is interactive, and requires a graphical session for its interactive R
plot for each sample and lineage.  (This is the `sonar_module_2_id_div_island`
rule in `workflow/sonar.smk`.)  If `--config sonar_member_list=shortcut` is
used, the manual step is skipped and sequences are selected from a prepared
[table of IDs](workflow/lineage_member_ids.csv).

Configuration options that can be set with `--config key=value`, with defaults
listed first:

 * `rundata`: `sra` to download raw read files from SRA, or a file path to
   process Illumina run directories from
 * `germline`: `germline-genbank` to use readymade IgDiscover outputs from
   GenBank, or `germline` to create them from scratch first
 * `sonar_member_list`: `shortcut` to extract lineage members from repertoire
   results with prepared lists from our own manually-guided identification
   (with confirmation that the sequences match), or `full` to run that
   manually-guided process here using SONAR's island selection procedure

[Snakemake]: https://snakemake.readthedocs.io
[Python]: https://www.python.org
[conda]: https://conda.io
[Singularity]: https://github.com/sylabs/singularity
[Docker]: https://www.docker.com/
[igseqhelper]: https://github.com/shawhahnlab/igseqhelper
[igseq]: https://github.com/shawhahnlab/igseq
[mamba]: https://mamba.readthedocs.io
[bcl2fastq]: https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html
[cutadapt]: https://github.com/marcelm/cutadapt
[PEAR]: https://cme.h-its.org/exelixis/web/software/pear/doc.html
[V2 apex antibody structure paper]: https://doi.org/10.1084/jem.20250638
