# nf-core/hic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.3.0dev - [date]

* New `--digestion` parameter to automatically set the restriction_site and ligation_site motifs
* New `--keep_multi` and `keep_dup` options. Default: false
* Template update for nf-core/tools v1.11
* Minor fix to summary log messages in pipeline header

### `Fixed`

* Fix recurrent bug in input file extension (#86)
* Fix bug in `--bin_size` parameter (#85)
* `min_mapq` is ignored if `--keep_multi` is used

### Deprecated

* `--rm_dup` and `rm_multi` are replaced by `--keep_dup` and `--keep_multi`

## v1.2.2 - 2020-09-02

### `Added`

* Template update for nf-core/tools v1.10.2
* Add the `--fastq_chunks_size` to specify the number of reads per chunks if split_fastq is true

### `Fixed`

* Bug in `--split_fastq` option not recognized

## v1.2.1 - 2020-07-06

### `Fixed`

* Fix issue with `--fasta` option and `.fa` extension (#66)

## v1.2.0 - 2020-06-18
