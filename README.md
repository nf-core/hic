# ![nf-core/hic](docs/images/nfcore-hic_logo.png)

**Analysis of Chromosome Conformation Capture data (Hi-C)**.

[![GitHub Actions CI Status](https://github.com/nf-core/hic/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/hic/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/hic/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/hic/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/hic.svg)](https://hub.docker.com/r/nfcore/hic)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2669513.svg)](https://doi.org/10.5281/zenodo.2669513)

## Introduction

This pipeline is based on the
[HiC-Pro workflow](https://github.com/nservant/HiC-Pro).
It was designed to process Hi-C data from raw FastQ files (paired-end Illumina
data) to normalized contact maps.
The current version supports most protocols, including digestion protocols as
well as protocols that do not require restriction enzymes such as DNase Hi-C.
In practice, this workflow was successfully applied to many data-sets including
dilution Hi-C, in situ Hi-C, DNase Hi-C, Micro-C, capture-C, capture Hi-C or
HiChip data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool
to run tasks across multiple compute infrastructures in a very portable manner.
It comes with docker / singularity containers making installation trivial and
results highly reproducible.

## Pipeline summary

1. Mapping using a two steps strategy to rescue reads spanning the ligation
sites (bowtie2)
2. Detection of valid interaction products
3. Duplicates removal
4. Create genome-wide contact maps at various resolution
5. Contact maps normalization using the ICE algorithm (iced)
6. Quality controls and report (MultiQC)
7. Addition export for visualisation and downstream analysis (cooler)

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/)
or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/)
for full pipeline reproducibility (please only use [`Conda`](https://conda.io/miniconda.html)
as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run nf-core/hic -profile test,<docker/singularity/conda/institute>
```

> Please check [nf-core/configs](https://github.com/nf-core/configs#documentation)
to see if a custom config file to run nf-core pipelines already exists for your Institute.
If so, you can simply use `-profile <institute>` in your command.
This will enable either `docker` or `singularity` and set the appropriate execution
settings for your local compute environment.

iv. Start running your own analysis!

```bash
nextflow run nf-core/hic -profile <docker/singularity/conda/institute> --reads '*_R{1,2}.fastq.gz' --genome GRCh37
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The nf-core/hic pipeline comes with documentation about the pipeline,
found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

For further information or help, don't hesitate to get in touch on
[Slack](https://nfcore.slack.com/channels/hic).
You can join with [this invite](https://nf-co.re/join/slack).

## Credits

nf-core/hic was originally written by Nicolas Servant.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on
[Slack](https://nfcore.slack.com/channels/hic) (you can join with
[this invite](https://nf-co.re/join/slack)).

## Citation

If you use nf-core/hic for your analysis, please cite it using the following
doi: [10.5281/zenodo.2669513](https://doi.org/10.5281/zenodo.2669513)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg,
Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13.
doi:[10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).  
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)
