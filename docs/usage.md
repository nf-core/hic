# nf-core/hic: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/hic/usage](https://nf-co.re/hic/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The `nf-core-hic` pipeline is designed to work only with paired-end data. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
SAMPLE_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
SAMPLE_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/hic --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile.
See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/hic -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/hic
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/hic releases page](https://github.com/nf-core/hic/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give
configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from
[https://github.com/nf-core/configs](https://github.com/nf-core/configs)
when it runs, making multiple config profiles for various institutional
clusters available at run time.
For more information and to see if your system is available in these
configs please see
the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` -
the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite
earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command).
See the [nf-core website documentation](https://nf-co.re/usage/configuration)
for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on
[Slack](https://nf-co.re/join/slack) on the
[`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs.
The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal
so that the workflow does not stop if you log out of your session. The logs are
saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached
session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted
your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a
large amount of memory.
We recommend adding the following line to your environment to limit this
(typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Use case

### Hi-C digestion protocol

Here is an command line example for standard DpnII digestion protocols.
Alignment will be performed on the `mm10` genome with default parameters.
Multi-hits will not be considered and duplicates will be removed.
Note that by default, no filters are applied on DNA and restriction fragment sizes.

```bash
nextflow run main.nf --input './*_R{1,2}.fastq.gz' --genome 'mm10' --digestion 'dnpii'
```

### DNase Hi-C / Micro-C protocol

Here is an command line example for DNase or Micro-C protocol.
Alignment will be performed on the `mm10` genome with default paramters.
Multi-hits will not be considered and duplicates will be removed.
Contacts involving fragments separated by less than 1000bp will be discarded.

```bash
nextflow run main.nf --input './*_R{1,2}.fastq.gz' --genome 'mm10' --no_digestion --min_cis 1000
```

## Inputs

### `--input`

Use this to specify the location of your input FastQ files. For example:

```bash
--input 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}`
   notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

Note that the Hi-C data analysis workflow requires paired-end data.

## Reference genomes

The pipeline config files come bundled with paths to the Illumina iGenomes reference
index files. If running with docker or AWS, the configuration is set up to use the
[AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)

There are many different species supported in the iGenomes references. To run
the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the
[iGenomes config file](https://github.com/nf-core/hic/blob/master/conf/igenomes.config).

### `--fasta`

If you prefer, you can specify the full path to your reference genome when you
run the pipeline:

```bash
--fasta '[path to Fasta reference]'
```

### `--bwt2_index`

The bowtie2 indexes are required to align the data with the HiC-Pro workflow. If the
`--bwt2_index` is not specified, the pipeline will either use the iGenomes
bowtie2 indexes (see `--genome` option) or build the indexes on-the-fly
(see `--fasta` option)

```bash
--bwt2_index '[path to bowtie2 index]'
```

### `--chromosome_size`

The Hi-C pipeline also requires a two-column text file with the
chromosome name and the chromosome size (tab-separated).
If not specified, this file will be automatically created by the pipeline.
In the latter case, the `--fasta` reference genome has to be specified.

```bash
   chr1    249250621
   chr2    243199373
   chr3    198022430
   chr4    191154276
   chr5    180915260
   chr6    171115067
   chr7    159138663
   chr8    146364022
   chr9    141213431
   chr10   135534747
   (...)
```

```bash
--chromosome_size '[path to chromosome size file]'
```

### `--restriction_fragments`

Finally, Hi-C experiments based on restriction enzyme digestion require a BED
file with coordinates of restriction fragments.

```bash
   chr1   0       16007   HIC_chr1_1    0   +
   chr1   16007   24571   HIC_chr1_2    0   +
   chr1   24571   27981   HIC_chr1_3    0   +
   chr1   27981   30429   HIC_chr1_4    0   +
   chr1   30429   32153   HIC_chr1_5    0   +
   chr1   32153   32774   HIC_chr1_6    0   +
   chr1   32774   37752   HIC_chr1_7    0   +
   chr1   37752   38369   HIC_chr1_8    0   +
   chr1   38369   38791   HIC_chr1_9    0   +
   chr1   38791   39255   HIC_chr1_10   0   +
   (...)
```

If not specified, this file will be automatically created by the pipeline.
In this case, the `--fasta` reference genome will be used.
Note that the `--digestion` or `--restriction_site` parameter is mandatory to create this file.

## Hi-C specific options

The following options are defined in the `nextflow.config` file, and can be
updated either using a custom configuration file (see `-c` option) or using
command line parameters.

### HiC-pro mapping

The reads mapping is currently based on the two-steps strategy implemented in
the HiC-pro pipeline. The idea is to first align reads from end-to-end.
Reads that do not align are then trimmed at the ligation site, and their 5'
end is re-aligned to the reference genome.
Note that the default options are quite stringent, and can be updated according
to the reads quality or the reference genome.

#### `--bwt2_opts_end2end`

Bowtie2 alignment option for end-to-end mapping.
Default: '--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end
--reorder'

```bash
--bwt2_opts_end2end '[Options for bowtie2 step1 mapping on full reads]'
```

#### `--bwt2_opts_trimmed`

Bowtie2 alignment option for trimmed reads mapping (step 2).
Default: '--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end
--reorder'

```bash
--bwt2_opts_trimmed '[Options for bowtie2 step2 mapping on trimmed reads]'
```

#### `--min_mapq`

Minimum mapping quality. Reads with lower quality are discarded. Default: 10

```bash
--min_mapq '[Minimum quality value]'
```

### Digestion Hi-C

#### `--digestion`

This parameter allows to automatically set the `--restriction_site` and
`--ligation_site` parameter according to the restriction enzyme you used.
Available keywords are 'hindiii', 'dpnii', 'mboi', 'arima'.

```bash
--digestion 'hindiii'
```

#### `--restriction_site`

If the restriction enzyme is not available through the `--digestion`
parameter, you can also define manually the restriction motif(s) for
Hi-C digestion protocol.
The restriction motif(s) is(are) used to generate the list of restriction fragments.
The precise cutting site of the restriction enzyme has to be specified using
the '^' character. Default: 'A^AGCTT'
Here are a few examples:

- MboI: ^GATC
- DpnII: ^GATC
- HindIII: A^AGCTT
- ARIMA kit: ^GATC,G^ANTC

Note that multiples restriction motifs can be provided (comma-separated) and
that 'N' base are supported.

```bash
--restriction_size '[Cutting motif]'
```

#### `--ligation_site`

Ligation motif after reads ligation. This motif is used for reads trimming and
depends on the fill in strategy.
Note that multiple ligation sites can be specified (comma-separated) and that
'N' base is interpreted and replaced by 'A','C','G','T'.
Default: 'AAGCTAGCTT'

```bash
--ligation_site '[Ligation motif]'
```

Exemple of the ARIMA kit: GATCGATC,GANTGATC,GANTANTC,GATCANTC

### DNAse/Micro-C

#### `--no_digestion`

In DNAse Hi-C mode, all options related to digestion Hi-C
(see previous section) are ignored.
In this case, it is highly recommended to use the `--min_cis_dist` parameter
to remove spurious ligation products.

```bash
--no_digestion
```

### HiC-pro processing

#### `--min_restriction_fragment_size`

Minimum size of restriction fragments to consider for the Hi-C processing.
Default: '0' - no filter

```bash
--min_restriction_fragment_size '[numeric]'
```

#### `--max_restriction_fragment_size`

Maximum size of restriction fragments to consider for the Hi-C processing.
Default: '0' - no filter

```bash
--max_restriction_fragment_size '[numeric]'
```

#### `--min_insert_size`

Minimum reads insert size. Shorter 3C products are discarded.
Default: '0' - no filter

```bash
--min_insert_size '[numeric]'
```

#### `--max_insert_size`

Maximum reads insert size. Longer 3C products are discarded.
Default: '0' - no filter

```bash
--max_insert_size '[numeric]'
```

#### `--min_cis_dist`

Filter short range contact below the specified distance.
Mainly useful for DNAse/Micro-C. Default: '0'

```bash
--min_cis_dist '[numeric]'
```

#### `--keep_dups`

If specified, duplicate reads are not discarded before building contact maps.

```bash
--keep_dups
```

#### `--keep_multi`

If specified, reads that aligned multiple times on the genome are not discarded.
Note the default mapping options are based on random hit assignment, meaning
that only one position is kept per read.
Note that in this case the `--min_mapq` parameter is ignored.

```bash
--keep_multi
```

## Genome-wide contact maps

Once the list of valid pairs is available, the standard is now to move on the `cooler`
framework to build the raw and balanced contact maps in txt and (m)cool formats.

### `--bin_size`

Resolution of contact maps to generate (comma-separated).
Default:'1000000,500000'

```bash
--bins_size '[string]'
```

### `--res_zoomify`

Define the maximum resolution to reach when zoomify the cool contact maps.
Default:'5000'

```bash
--res_zoomify '[string]'
```

### HiC-Pro contact maps

By default, the contact maps are now generated with the `cooler` framework.
However, for backward compatibility, the raw and normalized maps can still be generated
by HiC-pro if the `--hicpro_maps` parameter is set.

#### `--hicpro_maps`

If specified, the raw and ICE normalized contact maps will be generated by HiC-Pro.

```bash
--hicpro_maps
```

#### `--ice_max_iter`

Maximum number of iteration for ICE normalization.
Default: 100

```bash
--ice_max_iter '[numeric]'
```

#### `--ice_filer_low_count_perc`

Define which percentage of bins with low counts should be forced to zero.
Default: 0.02

```bash
--ice_filter_low_count_perc '[numeric]'
```

#### `--ice_filer_high_count_perc`

Define which percentage of bins with low counts should be discarded before
normalization. Default: 0

```bash
--ice_filter_high_count_perc '[numeric]'
```

#### `--ice_eps`

The relative increment in the results before declaring convergence for ICE
normalization. Default: 0.1

```bash
--ice_eps '[numeric]'
```

## Downstream analysis

### Additional quality controls

#### `--res_dist_decay`

Generates distance vs Hi-C counts plots at a given resolution using `HiCExplorer`.
Several resolutions can be specified (comma-separeted). Default: '250000'

```bash
--res_dist_decay '[string]'
```

### Compartment calling

Call open/close compartments for each chromosome, using the `cooltools` command.

#### `--res_compartments`

Resolution to call the chromosome compartments (comma-separated).
Default: '250000'

```bash
--res_compartments '[string]'
```

### TADs calling

#### `--tads_caller`

TADs calling can be performed using different approaches.
Currently available options are `insulation` and `hicexplorer`.
Note that all options can be specified (comma-separated).
Default: 'insulation'

```bash
--tads_caller '[string]'
```

#### `--res_tads`

Resolution to run the TADs calling analysis (comma-separated).
Default: '40000,20000'

```bash
--res_tads '[string]'
```

## Inputs/Outputs

### `--split_fastq`

By default, the nf-core Hi-C pipeline expects one read pairs per sample.
However, for large Hi-C data processing single fastq files can be very
time consuming.
The `--split_fastq` option allows to automatically split input read pairs
into chunks of reads of size `--fastq_chunks_size` (Default : 20000000). In this case, all chunks will be processed in parallel
and merged before generating the contact maps, thus leading to a significant
increase of processing performance.

```bash
--split_fastq --fastq_chunks_size '[numeric]'
```

### `--save_reference`

If specified, annotation files automatically generated from the `--fasta` file
are exported in the results folder. Default: false

```bash
--save_reference
```

### `--save_aligned_intermediates`

If specified, all intermediate mapping files are saved and exported in the
results folder. Default: false

```bash
--save_aligned_inermediates
```

### `--save_interaction_bam`

If specified, write a BAM file with all classified reads (valid pairs,
dangling end, self-circle, etc.) and its tags.

```bash
--save_interaction_bam
```

## Skip options

### `--skip_maps`

If defined, the workflow stops with the list of valid interactions, and the
genome-wide maps are not built. Useful for capture-C analysis. Default: false

```bash
--skip_maps
```

### `--skip_balancing`

If defined, the contact maps normalization is not run on the raw contact maps.
Default: false

```bash
--skip_balancing
```

### `--skip_cool`

If defined, cooler files are not generated. Default: false

```bash
--skip_cool
```

### `--skip_dist_decay`

Do not run distance decay plots. Default: false

```bash
--skip_dist_decay
```

### `--skip_compartments`

Do not call compartments. Default: false

```bash
--skip_compartments
```

### `--skip_tads`

Do not call TADs. Default: false

```bash
--skip_tads
```

### `--skip_multiQC`

If defined, the MultiQC report is not generated. Default: false

```bash
--skip_multiQC
```
