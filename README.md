# nf-core/hic

**Analysis of Chromosome Conformation Capture data (Hi-C)**.

[![Build Status](https://travis-ci.com/nf-core/hic.svg?branch=master)](https://travis-ci.com/nf-core/hic)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/hic.svg)](https://hub.docker.com/r/nfcore/hic)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
This pipeline is based on the [HiC-Pro workflow](https://github.com/nservant/HiC-Pro).
It was designed to process Hi-C data from raw fastq files (paired-end Illumina data) to normalized contact maps.
The current version supports most protocols, including digestion protocols as well as protocols that do not require restriction enzymes such as DNase Hi-C.
In practice, this workflow was successfully applied to many data-sets including dilution Hi-C, in situ Hi-C, DNase Hi-C, Micro-C, capture-C, capture Hi-C or HiChip data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.

### Pipeline summary
1. Mapping using a two steps strategy to rescue reads spanning the ligation sites (bowtie2)
2. Detection of valid interaction products
3. Duplicates removal
4. Create genome-wide contact maps at various resolution
5. Contact maps normalization using the ICE algorithm (iced)
6. Quality controls and report (MultiQC)
7. Addition export for visualisation and downstream analysis (cooler)

### Documentation
The nf-core/hic pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
    * [Reference genomes](docs/configuration/reference_genomes.md)  
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
nf-core/hic was originally written by Nicolas Servant.
