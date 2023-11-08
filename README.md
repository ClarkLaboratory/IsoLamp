# IsoLamp: isoform identification and quantification from long-read amplicon sequencing data

LOGO.

[![Install](https://img.shields.io/badge/Install-Github-brightgreen)](#installation)
[![Generic badge](https://img.shields.io/badge/Language-Bash-<COLOR>.svg)](https://shields.io/)
[![Generic badge](https://img.shields.io/badge/Publication-01.xx-<COLOR>.svg)](https://shields.io/)

**IsoLamp is a bash pipeline for the identification of known and novel isoforms from targeted amplicon long-read sequencing data generated with Oxford Nanopore Technologies.**

PIPELINE.

## Contents

- [Installation](#installation)
- [General Command Line Usage](#general-command-line-usage)
- [SIRV Test Dataset](#SIRV-test-dataset)
- [Dependencies](#dependencies)
- [Input Data](#input-data)
- [Parameters](#parameters)
- [Output](#Output)
- [Visualisation with IsoVis](#Visualisation-with-IsoVis)
- [Troubleshooting](#Troubleshooting)

## Installation
Download the pipeline with Git:
```
git clone https://github.com/youyupei/discoAnt.git
```

Run the setup script:
```
cd IsoLamp
bash scripts/IsoLamp_setup
```
**Note:** depending on your operating system you may need to edit the setup script to either 'conda activate' or 'source activate'.

## General Command Line Usage
We suggest adding IsoLamp to your PATH so it can be run from any directory:
```
echo 'export PATH=~/path/to/IsoLamp:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

Run IsoLamp (after editing the parameters file):
```
conda activate IsoLamp

IsoLamp parameters.ini 

# alternatively, if not in PATH:
./IsoLamp parameters.ini # if running from same directory
~/path/to/IsoLamp/IsoLamp parameters.ini # provide full path when running from other directories
```

## SIRV Test Dataset
Test the installation on the provided SIRV data:
```
conda activate IsoLamp
IsoLamp sirv_test_data/sirv_params.ini
```
This should produce a folder called 'SIRV5_test' which contains the expected output of the pipeline.

## Dependencies
All dependencies and R libraries are automatically built in a conda environment when the setup script is executed.

Packages:
  - python=3.7
  - pandas
  - pysam
  - numpy
  - tqdm
  - BBMap
  - bedtools
  - samtools
  - salmon=0.14.2
  - gffread
  - gffcompare
  - minimap2


  R and libraries:
  - R>=4.3
  - Bioconductor (R package)
  - bambu>=3.2.4 (Bioconductor)
  - optparse
  - dplyr
  - factominer
  - factoextra
  - ggrepel
  - ggplot2
  - reshape
  - rstatix
  - purrr

## Input Data
The pipeline is designed to run on barcoded data produced by Oxford Nanopore Technologies, typically from amplicon sequencing of a gene of interest. This is typically a directory containing subdirectories titled 'barcode01', 'barcode02' etc. Each barcode directory contains either single or multiple FASTA or FASTQ files of reads. If you have multiple experiments amplifying different genes, run the pipeline once for each gene.

Please see the sirv_test_data/test_fasta directory for an example of the structure required. If you don't have barcoded data, please create a subdirectory with your single sample name and place all read FASTA or FASTQ files inside this directory (nanopore_data/barcode01/reads.fastq, nanopore_data/barcode02/reads.fastq, ...).

**Note:** the 'barcode' name is arbitrary, your subdirectories can be renamed to sample names or any identifiers.

## Parameters
Basic parameters file:
```
OUTPUT_NAME=CLCN3
ENSG_ID=ENSG00000109572
READS="path/to/sample/fa_folders"
GENOME="path/to/hg38.fa"
ANNOTATION="path/to/gencode.v44.annotation.gtf"
```
Full parameters file (with comments):
```
## Experiment Info ##
# a folder will be created with OUTPUT_NAME and this is used as a prefix for files
OUTPUT_NAME=CLCN3
# the ENSEMBL ID is used to filter reference files to increase speed
ENSG_ID=ENSG00000109572


## Paths to data ##
# standard output for barcoded amplicon sequencing is multiple barcode/sample folers containing fastq or fasta files
# supports .fasta or .fastq
READS="path/to/sample/fa_folders"
# path to reference genome 
GENOME="path/to/hg38.fa"
# path to reference annotation
ANNOTATION="path/to/gencode.v44.annotation.gtf"

# path to a CSV of sample IDs and group IDs
grouping_data="path/to/group_ids.csv" # default is NULL, leave blank for default

## Minimum read count per isoform threshold ##
read_count_minimum=5 # default is 5, leave blank for default
samples_minimum=3 # default is half the total number of samples (rounded down to integer), leave blank for default

## Downsampling to be set to TRUE or FALSE ##
# Downsampling results in the same number of reads for every sample, reducing inter-sample bias
downsampling=TRUE # default is TRUE
number_reads_downsample=8000 # default is 8000, ignored if downsampling=FALSE

## Primer site filter to be set to TRUE or FALSE ##
# Recommended to increase accuracy and reduce isoforms to only those amplified with primers
primer_site_based_filter=TRUE # default is FALSE
# Primer coordinates for forward and reverse, as bed files
# Only checked if primer filtering is TRUE
forward_primers="path/to/forward.bed" # default is NULL
reverse_primers="path/to/reverse.bed" # default is NULL

## Extract high quality splice junctions to be set to TRUE or FALSE ##
# Recommended to increase novel isoform accuracy
extract_high_quality_SJs=TRUE # default is FALSE
# minimum junction alignment quality
JAQ=0.8
# window (nt) upstream and downstream of splice junction
junction_window=25 

## Other options ##
# minimap2 -G intron length, leave blank for default
max_intron_length=400 
# bambu parameters, leave blank for default
bambu_ndr=1 
bambu_min_gene_fraction=0.001
```

Detailed descriptions of parameters:
|Parameter|Type|Default|Description| 
|---|---|---|---|
|OUTPUT_NAME|required|NA|A directory will be created with this name, we suggest using the gene name.|
|ENSG_ID|required|NA|The ENSEMBL gene ID. This is used to filter the reference genome and annotation. An ENSEMBL ID is also required for visualisation with IsoVis. If this ID not available, use the gene name as it appears in the GTF (as in the 'sirv_test_data/sirv_params.ini' file).|
|READS|required|NA|Path to the top level directory of sample/barcode directories.|
|GENOME|required|NA|The reference genome FASTA.|
|ANNOTATION|required|NA|The reference annotation GTF.|
|grouping_data|optional|NULL|A CSV file of sample and group names used for plotting. The first row must be 'sample,group'. See example of this below and in 'sirv_test_data/sirv_grouping.csv'|
|read_count_minimum|optional|5|Known and novel isoforms must meet this threshold in at least 'samples_minimum' number of samples.|
|samples_minimum|optional|NA|Known and novel isoforms must meet the read count minimum in this many samples. Default value caclulated as number of samples/2 rounded down to nearest integer.|
|downsampling|optional|TRUE|Whether each sample/barcode should be downsampled to a consistent number of reads.|
|number_reads_downsample|optional|8000|Number of reads to downsample each sample/barcode.|
|primer_site_based_filter|optional|FALSE|Whether to remove isoforms that do not overlap the primers used to perform amplicon sequencing. We recommend using this option.|
|forward_primers|optional|NULL|BED file of forward primers. Only checked if primer_site_based_filter is TRUE.|
|reverse_primers|optional|NULL|BED file of reverse primers. Only checked if primer_site_based_filter is TRUE.|
|extract_high_quality_SJs|optional|FALSE|Whether to extract reads with high quality splice junctions, and run Bambu on these reads. We recommend using this option to increase accuracy in calling novel isoforms.|
|JAQ|optional|0.8|The minimum junction alignment quality every junction must meet in order for a read to be included. Only checked if extract_high_quality_SJs is TRUE.|
|junction_window|optional|25|The nt distance upstream and downstream of splice junctions to calculate the JAQ. Only checked if extract_high_quality_SJs is TRUE.|
|max_intron_length|optional|400|Controls the '-G' flag in minimap2 as some complex genes have long introns.|
|bambu_ndr|optional|1|Controls the bambu 'NDR' option. We set this to 1 by default to return all possible novel isoforms and filter these downstream.|
|bambu_min_gene_fraction|optional|0.001|Controls the bambu 'min.readFractionByGene' option.|

Example grouping_data file:
```
sample,group
barcode01,group1
barcode02,group1
barcode03,group2
barcode04,group2
```

## Output
The main output of IsoLamp includes:
  - a basic text report
  - annotation of known and novel isoforms as a GTF
  - information on novel isoforms classes from Gffcompare and Bambu
  - quantification of known and novel isoforms (counts and propoprtions) as CSVs
  - a PCA plot of samples/barcodes
  - an accuracy plot of all input reads
  - if applicable, a CSV output of t-test results comparing isoform proportions between groups

## Visualisation with IsoVis
The results from the IsoLamp pipeline can be visualised using IsoVis: https://isomix.org/isovis.

The isoform annotation GTF (and optionally isoform counts CSV) can be directly uploaded to IsoVIs.

## Troubleshooting






