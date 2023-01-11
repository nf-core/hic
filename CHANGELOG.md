# nf-core/hic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - 2023-01-12

### `Added`

- DSL2 version of nf-core-hic pipeline
- Add full test dataset (#80)
- Replace local modules by the cooler nf-core module

### `Fixed`

- Fix error in the Arima preset (#127)

## v1.3.1 - 2021-09-25

### `Fixed`

- Fix bug in conda environment for cooltools (#109)

## v1.3.0 - 2021-05-22

- Change the `/tmp/` folder to `./tmp/` folder so that all tmp files are now in the work directory (#24)
- Add `--hicpro_maps` options to generate the raw and normalized HiC-Pro maps. The default is now to use cooler
- Add chromosome compartments calling with cooltools (#53)
- Add HiCExplorer distance decay quality control (#54)
- Add HiCExplorer TADs calling (#55)
- Add insulation score TADs calling (#55)
- Generate cooler/txt contact maps
- Normalize Hi-C data with cooler instead of iced
- New `--digestion` parameter to automatically set the restriction_site and ligation_site motifs
- New `--keep_multi` and `keep_dup` options. Default: false
- Template update for nf-core/tools
- Minor fix to summary log messages in pipeline header

### `Fixed`

- Fix bug in stats report which were not all correcly exported in the results folder
- Fix recurrent bug in input file extension (#86)
- Fix bug in `--bin_size` parameter (#85)
- `--min_mapq` is ignored if `--keep_multi` is used

### `Deprecated`

- `--rm_dup` and `--rm_multi` are replaced by `--keep_dups` and `--keep_multi`

## v1.2.2 - 2020-09-02

### `Added`

- Template update for nf-core/tools v1.10.2
- Add the `--fastq_chunks_size` to specify the number of reads per chunks if split_fastq is true

### `Fixed`

- Bug in `--split_fastq` option not recognized

## v1.2.1 - 2020-07-06

### `Fixed`

- Fix issue with `--fasta` option and `.fa` extension (#66)

## v1.2.0 - 2020-06-18

### `Added`

- Bump v1.2.0
- Merge template nf-core 1.9
- Move some options to camel_case
- Update python scripts for python3
- Update conda environment file
  - python base `2.7.15` > `3.7.6`
  - pip `19.1` > `20.0.1`
  - scipy `1.2.1` > `1.4.1`
  - numpy `1.16.3` > `1.18.1`
  - bx-python `0.8.2` > `0.8.8`
  - pysam `0.15.2` > `0.15.4`
  - cooler `0.8.5` > `0.8.6`
  - multiqc `1.7` > `1.8`
  - iced `0.5.1` > `0.5.6`
  - _*New*_ pymdown-extensions `7.1`
  - _*New*_ hicexplorer `3.4.3`
  - _*New*_ bioconductor-hitc `1.32.0`
  - _*New*_ r-optparse `1.6.6`
  - _*New*_ ucsc-bedgraphtobigwig `377`
  - _*New*_ cython `0.29.19`
  - _*New*_ cooltools `0.3.2`
  - _*New*_ fanc `0.8.30`
  - _*Removed*_ r-markdown

### `Fixed`

- Fix error in doc for Arima kit usage
- Sort output of `get_valid_interaction` process as the input files of `remove_duplicates`
  are expected to be sorted (sort -m)

### `Deprecated`

- Command line options converted to `camel_case`:
  - `--skipMaps` > `--skip_maps`
  - `--skipIce` > `--skip_ice`
  - `--skipCool` > `--skip_cool`
  - `--skipMultiQC` > `--skip_multiqc`
  - `--saveReference` > `--save_reference`
  - `--saveAlignedIntermediates` > `--save_aligned_intermediates`
  - `--saveInteractionBAM` > `--save_interaction_bam`

## v1.1.1 - 2020-04-02

### `Fixed`

- Fix bug in tag. Remove '['

## v1.1.0 - 2019-10-15

### `Added`

- Update hicpro2higlass with `-p` parameter
- Support 'N' base motif in restriction/ligation sites
- Support multiple restriction enzymes/ligattion sites (comma separated) ([#31](https://github.com/nf-core/hic/issues/31))
- Add --saveInteractionBAM option
- Add DOI ([#29](https://github.com/nf-core/hic/issues/29))
- Update manual ([#28](https://github.com/nf-core/hic/issues/28))

### `Fixed`

- Fix bug for reads extension `_1`/`_2` ([#30](https://github.com/nf-core/hic/issues/30))

## v1.0 - [2019-05-06]

Initial release of nf-core/hic, created with the [nf-core](http://nf-co.re/) template.

### `Added`

First version of nf-core Hi-C pipeline which is a Nextflow implementation of
the [HiC-Pro pipeline](https://github.com/nservant/HiC-Pro/).
Note that all HiC-Pro functionalities are not yet all implemented.
The current version supports most protocols including Hi-C, in situ Hi-C,
DNase Hi-C, Micro-C, capture-C or HiChip data.

In summary, this version allows :

- Automatic detection and generation of annotation files based on igenomes
  if not provided.
- Two-steps alignment of raw sequencing reads
- Reads filtering and detection of valid interaction products
- Generation of raw contact matrices for a set of resolutions
- Normalization of the contact maps using the ICE algorithm
- Generation of cooler file for visualization on [higlass](https://higlass.io/)
- Quality report based on HiC-Pro MultiQC module
