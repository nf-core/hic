# nf-core/hic: Changelog

## v1.0dev - 2019-04-09

First version of nf-core-hic pipeline which is a Nextflow implementation of the [HiC-Pro pipeline](https://github.com/nservant/HiC-Pro/).
Note that all HiC-Pro functionalities are not yet all implemented. The current version is designed for protocols based on restriction enzyme digestion.

In summary, this version allows :

* Automatic detection and generation of annotation files based on igenomes if not provided.
* Two-steps alignment of raw sequencing reads
* Reads filtering and detection of valid interaction products
* Generation of raw contact matrices for a set of resolutions
* Normalization of the contact maps using the ICE algorithm
* Generation of cooler file for visualization on [higlass](https://higlass.io/)
* Quality report based on HiC-Pro MultiQC module
