# ![nf-core/hic](docs/images/nfcore-hic_logo.png)

**Analysis of Chromosome Conformation Capture data (Hi-C)**.

[![Build Status](https://travis-ci.com/nf-core/hic.svg?branch=master)](https://travis-ci.com/nf-core/hic)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/hic.svg)](https://hub.docker.com/r/nfcore/hic)
![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2669513.svg)](https://doi.org/10.5281/zenodo.2669513)

## Introduction

This pipeline is based on the
[HiC-Pro workflow](https://github.com/nservant/HiC-Pro).
It was designed to process Hi-C data from raw fastq files (paired-end Illumina
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

## Documentation

The nf-core/hic pipeline comes with documentation about the pipeline, found in
the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

nf-core/hic was originally written by Nicolas Servant.

If you use nf-core/hic for your analysis, please cite it using the following
doi: [10.5281/zenodo.2669513](https://doi.org/10.5281/zenodo.2669513)

You can cite the `nf-core` pre-print as follows:
Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di
Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**.
*bioRxiv*. 2019. p. 610741.
[doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).
