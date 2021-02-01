# nf-core/hic: Output

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/hic/output](https://nf-co.re/hic/output)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.
The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [Reads alignment](#reads-alignment)
* [Valid pairs detection](#valid-pairs-detection)
* [Duplicates removal](#duplicates-removal)
* [Contact maps](#contact-maps)
* [MultiQC](#multiqc) - aggregate report and quality controls, describing
results of the whole pipeline
* [Export](#exprot) - additionnal export for compatibility with downstream
analysis tool and visualization

The current version is mainly based on the
[HiC-Pro](https://github.com/nservant/HiC-Pro) pipeline.
For details about the workflow, see
[Servant et al. 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0831-x)

## Reads alignment

Using Hi-C data, each reads mate has to be independantly aligned on the
reference genome.
The current workflow implements a two steps mapping strategy. First, the reads
are aligned using an end-to-end aligner.
Second, reads spanning the ligation junction are trimmmed from their 3' end,
and aligned back on the genome.
Aligned reads for both fragment mates are then paired in a single paired-end
BAM file.
Singletons are discarded, and multi-hits are filtered according to the
configuration parameters (`--rm-multi`).
Note that if the `--dnase` mode is activated, HiC-Pro will skip the second
mapping step.

**Output directory: `results/mapping`**

* `*bwt2pairs.bam` - final BAM file with aligned paired data
* `*.pairstat` - mapping statistics

if `--saveAlignedIntermediates` is specified, additional mapping file results
are available ;

* `*.bam` - Aligned reads (R1 and R2) from end-to-end alignment
* `*_unmap.fastq` - Unmapped reads after end-to-end alignment
* `*_trimmed.fastq` - Trimmed reads after end-to-end alignment
* `*_trimmed.bam` - Alignment of trimmed reads
* `*bwt2merged.bam` - merged BAM file after the two-steps alignment
* `*.mapstat` - mapping statistics per read mate

Usually, a high fraction of reads is expected to be aligned on the genome
(80-90%). Among them, we usually observed a few percent (around 10%) of step 2
aligned reads. Those reads are chimeric fragments for which we detect a
ligation junction. An abnormal level of chimeric reads can reflect a ligation
issue during the library preparation.
The fraction of singleton or multi-hits depends on the genome complexity and
the fraction of unmapped reads. The fraction of singleton is usually close to
the sum of unmapped R1 and R2 reads, as it is unlikely that both mates from the
same pair were unmapped.

## Valid pairs detection

Each aligned reads can be assigned to one restriction fragment according to the
reference genome and the digestion protocol.

Invalid pairs are classified as follow:

* Dangling end, i.e. unligated fragments (both reads mapped on the same
restriction fragment)
* Self circles, i.e. fragments ligated on themselves (both reads mapped on the
same restriction fragment in inverted orientation)
* Religation, i.e. ligation of juxtaposed fragments
* Filtered pairs, i.e. any pairs that do not match the filtering criteria on
inserts size, restriction fragments size
* Dumped pairs, i.e. any pairs for which we were not able to reconstruct the
ligation product.

Only valid pairs involving two different restriction fragments are used to
build the contact maps.
Duplicated valid pairs associated to PCR artefacts are discarded
(see `--rm_dup`).

In case of Hi-C protocols that do not require a restriction enzyme such as
DNase Hi-C or micro Hi-C, the assignment to a restriction is not possible
(see `--dnase`).
Short range interactions that are likely to be spurious ligation products
can thus be discarded using the `--min_cis_dist` parameter.

* `*.validPairs` - List of valid ligation products
* `*.DEpairs` - List of dangling-end products
* `*.SCPairs` - List of self-circle products
* `*.REPairs` - List of religation products
* `*.FiltPairs` - List of filtered pairs
* `*RSstat` - Statitics of number of read pairs falling in each category

The validPairs are stored using a simple tab-delimited text format ;

```bash
read name / chr_reads1 / pos_reads1 / strand_reads1 / chr_reads2 / pos_reads2 /
strand_reads2 / fragment_size / res frag name R1 / res frag R2 / mapping qual R1
/ mapping qual R2 [/ allele_specific_tag]
```

The ligation efficiency can be assessed using the filtering of valid and
invalid pairs. As the ligation is a random process, 25% of each valid ligation
class is expected. In the same way, a high level of dangling-end or self-circle
read pairs is associated with a low quality experiment, and reveals a problem
during the digestion, fill-in or ligation steps.

In the context of Hi-C protocol without restriction enzyme, this analysis step
is skipped. The aligned pairs are therefore directly used to generate the
contact maps. A filter of the short range contact (typically <1kb) is
recommanded as this pairs are likely to be self ligation products.

## Duplicates removal

Note that validPairs file are generated per reads chunck.
These files are then merged in the allValidPairs file, and duplicates are
removed if the `--rm_dup` parameter is used.

* `*allValidPairs` - combined valid pairs from all read chunks
* `*mergestat` - statistics about duplicates removal and valid pairs information

Additional quality controls such as fragment size distribution can be extracted
from the list of valid interaction products.
We usually expect to see a distribution centered around 300 pb which correspond
to the paired-end insert size commonly used.
The fraction of dplicates is also presented. A high level of duplication
indicates a poor molecular complexity and a potential PCR bias.
Finaly, an important metric is to look at the fraction of intra and
inter-chromosomal interactions, as well as long range (>20kb) versus short
range (<20kb) intra-chromosomal interactions.

## Contact maps

Intra et inter-chromosomal contact maps are build for all specified resolutions.
The genome is splitted into bins of equal size. Each valid interaction is
associated with the genomic bins to generate the raw maps.
In addition, Hi-C data can contain several sources of biases which has to be
corrected.
The current workflow uses the [ìced](https://github.com/hiclib/iced) and
[Varoquaux and Servant, 2018](http://joss.theoj.org/papers/10.21105/joss.01286)
python package which proposes a fast implementation of the original ICE
normalization algorithm (Imakaev et al. 2012), making the assumption of equal
visibility of each fragment.

* `*.matrix` - genome-wide contact maps
* `*_iced.matrix` - genome-wide iced contact maps

The contact maps are generated for all specified resolution
(see `--bin_size` argument)
A contact map is defined by :

* A list of genomic intervals related to the specified resolution (BED format).
* A matrix, stored as standard triplet sparse format (i.e. list format).

Based on the observation that a contact map is symmetric and usually sparse,
only non-zero values are stored for half of the matrix. The user can specified
if the 'upper', 'lower' or 'complete' matrix has to be stored. The 'asis'
option allows to store the contacts as they are observed from the valid pairs
files.

```bash
   A   B   10
   A   C   23
   B   C   24
   (...)
```

This format is memory efficient, and is compatible with several software for
downstream analysis.

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single
HTML report summarising all samples in your project. Most of the pipeline QC
results are visualised in the report and further statistics are available in
within the report data directory.

The pipeline has special steps which allow the software versions used to be
reported in the MultiQC output for future traceability.

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`,
  `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`,
  `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format:
  `results_description.html`.
