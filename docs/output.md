# nf-core/hic: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.
The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [From raw data to valid pairs](#from-raw-data-to-valid-pairs)
  - [HiC-Pro](#hicpro)
    - [Reads alignment](#reads-alignment)
    - [Valid pairs detection](#valid-pairs-detection)
    - [Duplicates removal](#duplicates-removal)
    - [Contact maps](#hicpro-contact-maps)
- [Hi-C contact maps](#hic-contact-maps)
- [Downstream analysis](#downstream-analysis)
  - [Distance decay](#distance-decay)
  - [Compartments calling](#compartments-calling)
  - [TADs calling](#tads-calling)
- [MultiQC](#multiqc) - aggregate report and quality controls, describing
  results of the whole pipeline
- [Export](#exprot) - additionnal export for compatibility with downstream
  analysis tool and visualization

## From raw data to valid pairs

### HiC-Pro

The current version is mainly based on the
[HiC-Pro](https://github.com/nservant/HiC-Pro) pipeline.
For details about the workflow, see
[Servant et al. 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0831-x)

#### Reads alignment

Using Hi-C data, each reads mate has to be independantly aligned on the
reference genome.
The current workflow implements a two steps mapping strategy. First, the reads
are aligned using an end-to-end aligner.
Second, reads spanning the ligation junction are trimmmed from their 3' end,
and aligned back on the genome.
Aligned reads for both fragment mates are then paired in a single paired-end
BAM file.
Singletons and low quality mapped reads are filtered (`--min_mapq`).
Note that if the `--dnase` mode is activated, HiC-Pro will skip the second
mapping step.

**Output directory: `results/hicpro/mapping`**

- `*bwt2pairs.bam` - final BAM file with aligned paired data

if `--save_aligned_intermediates` is specified, additional mapping file results
are available ;

- `*.bam` - Aligned reads (R1 and R2) from end-to-end alignment
- `*_unmap.fastq` - Unmapped reads after end-to-end alignment
- `*_trimmed.fastq` - Trimmed reads after end-to-end alignment
- `*_trimmed.bam` - Alignment of trimmed reads
- `*bwt2merged.bam` - merged BAM file after the two-steps alignment
- `*.mapstat` - mapping statistics per read mate

Usually, a high fraction of reads is expected to be aligned on the genome
(80-90%). Among them, we usually observed a few percent (around 10%) of step 2
aligned reads. Those reads are chimeric fragments for which we detect a
ligation junction. An abnormal level of chimeric reads can reflect a ligation
issue during the library preparation.
The fraction of singleton or low quality reads depends on the genome complexity and
the fraction of unmapped reads. The fraction of singleton is usually close to
the sum of unmapped R1 and R2 reads, as it is unlikely that both mates from the
same pair were unmapped.

#### Valid pairs detection with HiC-Pro

Each aligned reads can be assigned to one restriction fragment according to the
reference genome and the digestion protocol.

Invalid pairs are classified as follow:

- Dangling end, i.e. unligated fragments (both reads mapped on the same
  restriction fragment)
- Self circles, i.e. fragments ligated on themselves (both reads mapped on the
  same restriction fragment in inverted orientation)
- Religation, i.e. ligation of juxtaposed fragments
- Filtered pairs, i.e. any pairs that do not match the filtering criteria on
  inserts size, restriction fragments size
- Dumped pairs, i.e. any pairs for which we were not able to reconstruct the
  ligation product.

Only valid pairs involving two different restriction fragments are used to
build the contact maps.
Duplicated valid pairs associated to PCR artefacts are discarded
(see `--keep_dup` to not discard them).

In case of Hi-C protocols that do not require a restriction enzyme such as
DNase Hi-C or micro Hi-C, the assignment to a restriction is not possible
(see `--dnase`).
Short range interactions that are likely to be spurious ligation products
can thus be discarded using the `--min_cis_dist` parameter.

**Output directory: `results/hicpro/valid_pairs`**

- `*.validPairs` - List of valid ligation products
- `*.DEpairs` - List of dangling-end products
- `*.SCPairs` - List of self-circle products
- `*.REPairs` - List of religation products
- `*.FiltPairs` - List of filtered pairs
- `*RSstat` - Statitics of number of read pairs falling in each category

Of note, these results are saved only if `--save_pairs_intermediates` is used.  
The `validPairs` are stored using a simple tab-delimited text format ;

```bash
read name / chr_reads1 / pos_reads1 / strand_reads1 / chr_reads2 / pos_reads2 /
strand_reads2 / fragment_size / res frag name R1 / res frag R2 / mapping qual R1
/ mapping qual R2
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

#### Duplicates removal

Note that `validPairs` file are generated per reads chunck (and saved only if
`--save_pairs_intermediates` is specified).
These files are then merged in the `allValidPairs` file, and duplicates are
removed (see `--keep_dups` to disable duplicates filtering).

**Output directory: `results/hicpro/valid_pairs`**

- `*allValidPairs` - combined valid pairs from all read chunks

Additional quality controls such as fragment size distribution can be extracted
from the list of valid interaction products.
We usually expect to see a distribution centered around 300 pb which correspond
to the paired-end insert size commonly used.
The fraction of dplicates is also presented. A high level of duplication
indicates a poor molecular complexity and a potential PCR bias.
Finally, an important metric is to look at the fraction of intra and
inter-chromosomal interactions, as well as long range (>20kb) versus short
range (<20kb) intra-chromosomal interactions.

#### Pairs file

`.pairs` is a standard tabular format proposed by the 4DN Consortium
for storing DNA contacts detected in a Hi-C experiment
(see https://pairtools.readthedocs.io/en/latest/formats.html).
This format is the entry point of the downstream steps of the pipeline after
detection of valid pairs.

**Output directory: `results/hicpro/valid_pairs/pairix`**

- `*pairix` - compressed and indexed pairs file

#### Statistics

Various statistics files are generated all along the data processing.
All results are available in `results/hicpro/stats`.

**Output directory: `results/hicpro/stats`**

- \*mapstat - mapping statistics per read mate
- \*pairstat - R1/R2 pairing statistics
- \*RSstat - Statitics of number of read pairs falling in each category
- \*mergestat - statistics about duplicates removal and valid pairs information

#### Contact maps

Intra et inter-chromosomal contact maps are build for all specified resolutions.
The genome is splitted into bins of equal size. Each valid interaction is
associated with the genomic bins to generate the raw maps.
In addition, Hi-C data can contain several sources of biases which has to be
corrected.
The HiC-Pro workflow uses the [ìced](https://github.com/hiclib/iced) and
[Varoquaux and Servant, 2018](http://joss.theoj.org/papers/10.21105/joss.01286)
python package which proposes a fast implementation of the original ICE
normalization algorithm (Imakaev et al. 2012), making the assumption of equal
visibility of each fragment.

Importantly, the HiC-Pro maps are generated only if the `--hicpro_maps` option
is specified on the command line.

**Output directory: `results/hicpro/matrix`**

- `*.matrix` - genome-wide contact maps
- `*_iced.matrix` - genome-wide iced contact maps

The contact maps are generated for all specified resolutions
(see `--bin_size` argument).  
A contact map is defined by :

- A list of genomic intervals related to the specified resolution (BED format).
- A matrix, stored as standard triplet sparse format (i.e. list format).

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

## Hi-C contact maps

Contact maps are usually stored as simple txt (`HiC-Pro`), .hic (`Juicer/Juicebox`) and .(m)cool (`cooler/Higlass`) formats.
The .cool and .hic format are compressed and indexed and usually much more efficient that the txt format.  
In the current workflow, we propose to use the `cooler` format as a standard to build the raw and normalized maps
after valid pairs detection as it is used by several downstream analysis and visualization tools.

Raw contact maps are therefore in **`results/contact_maps/raw`** which contains the different maps in `txt` and `cool` formats, at various resolutions.
Normalized contact maps are stored in **`results/contact_maps/norm`** which contains the different maps in `txt`, `cool`, and `mcool` format.
The bin coordinates used for all resolutions are available in **`results/contact_maps/bins`**.

Note that `txt` contact maps generated with `cooler` are identical to those generated by `HiC-Pro`.
However, differences can be observed on the normalized contact maps as the balancing algorithm is not exactly the same.

## Downstream analysis

Downstream analysis are performed from `cool` files at specified resolution.

### Distance decay

The distance decay plot shows the relationship between contact frequencies and genomic distance. It gives a good indication of the compaction of the genome.
According to the organism, the slope of the curve should fit the expectation of polymer physics models.

The results generated with the `HiCExplorer hicPlotDistVsCounts` tool (plot and table) are available in the **`results/dist_decay/`** folder.

### Compartments calling

Compartments calling is one of the most common analysis which aims at detecting A (open, active) / B (close, inactive) compartments.
In the first studies on the subject, the compartments were called at high/medium resolution (1000000 to 250000) which is enough to call A/B comparments.
Analysis at higher resolution has shown that these two main types of compartments can be further divided into compartments subtypes.

Although different methods have been proposed for compartment calling, the standard remains the eigen vector decomposition from the normalized correlation maps.
Here, we use the implementation available in the [`cooltools`](https://cooltools.readthedocs.io/en/lates) package.

Results are available in **`results/compartments/`** folder and includes :

- `*cis.vecs.tsv`: eigenvectors decomposition along the genome
- `*cis.lam.txt`: eigenvalues associated with the eigenvectors

### TADs calling

TADs has been described as functional units of the genome.
While contacts between genes and regulatority elements can occur within a single TADs, contacts between TADs are much less frequent, mainly due to the presence of insulation protein (such as CTCF) at their boundaries. Looking at Hi-C maps, TADs look like triangles around the diagonal. According to the contact map resolutions, TADs appear as hierarchical structures with a median size around 1Mb (in mammals), as well as smaller structures usually called sub-TADs of smaller size.

TADs calling remains a challenging task, and even if many methods have been proposed in the last decade, little overlap have been found between their results.

Currently, the pipeline proposes two approaches :

- Insulation score using the [`cooltools`](https://cooltools.readthedocs.io/en/latest/cli.html#cooltools-diamond-insulation) package. Results are availabe in **`results/tads/insulation`**.
- [`HiCExplorer TADs calling`](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html). Results are available at **`results/tads/hicexplorer`**.

Usually, TADs results are presented as simple BED files, or bigWig files, with the position of boundaries along the genome.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
