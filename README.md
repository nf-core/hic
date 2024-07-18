<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-hic_logo_dark.png">
    <img alt="nf-core/hic" src="docs/images/nf-core-hic_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/hic/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/hic/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/hic/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/hic/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/hic/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.2669512-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.2669512)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/hic)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23hic-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/hic)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/hic** is a bioinformatics best-practice analysis pipeline for Analysis of Chromosome Conformation Capture data (Hi-C).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/hic/results).

## Pipeline summary

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Hi-C data processing
   1. [`HiC-Pro`](https://github.com/nservant/HiC-Pro)
      1. Mapping using a two steps strategy to rescue reads spanning the ligation
         sites ([`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
      2. Detection of valid interaction products
      3. Duplicates removal
      4. Generate raw and normalized contact maps ([`iced`](https://github.com/hiclib/iced))
      5. Generate `pairs` files for downstream analysis
   2. [`Pairtools`](https://github.com/open2c/pairtools)
      1. Mapping using [`BWA-mem`](https://github.com/lh3/bwa)
      2. Detection of valid interaction products with [`pairtools`](https://github.com/open2c/pairtools)
      3. Duplicates removal
      4. Generate `pairs` files for downstream analysis
3. Create genome-wide contact maps at various resolutions ([`cooler`](https://github.com/open2c/cooler))
4. Contact maps normalization using balancing algorithm ([`cooler`](https://github.com/open2c/cooler))
5. Export to various contact maps formats ([`HiC-Pro`](https://github.com/nservant/HiC-Pro), [`cooler`](https://github.com/open2c/cooler))
6. Quality controls ([`HiC-Pro`](https://github.com/nservant/HiC-Pro), [`HiCExplorer`](https://github.com/deeptools/HiCExplorer))
7. Compartments calling ([`cooltools`](https://cooltools.readthedocs.io/en/latest/), [`Calder2`](https://github.com/CSOgroup/CALDER2))
8. TADs calling ([`HiCExplorer`](https://github.com/deeptools/HiCExplorer), [`cooltools`](https://cooltools.readthedocs.io/en/latest/))
9. Quality control report ([`MultiQC`](https://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
HIC_ES_4,SRR5339783_1.fastq.gz,SRR5339783_2.fastq.gz
```

Each row represents a pair of fastq files (paired end).
Now, you can run the pipeline using:

```bash
nextflow run nf-core/hic \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --genome GRCh37 \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/hic/usage) and the [parameter documentation](https://nf-co.re/hic/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/hic/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/hic/output).

## Credits

nf-core/hic was originally written by Nicolas Servant.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#hic` channel](https://nfcore.slack.com/channels/hic) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/hic for your analysis, please cite it using the following doi: [10.5281/zenodo.2669512](https://doi.org/10.5281/zenodo.2669512)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
