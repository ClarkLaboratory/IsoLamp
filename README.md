# IsoLamp: isoform identification and quantification from long-read amplicon sequencing data

[![Install](https://img.shields.io/badge/Install-Github-brightgreen)](#installation)
[![Generic badge](https://img.shields.io/badge/Language-Bash-<COLOR>.svg)](https://shields.io/)
[![Generic badge](https://img.shields.io/badge/Publication-10.1101/2024.02.22.24303189-<COLOR>.svg)](https://www.medrxiv.org/content/10.1101/2024.02.22.24303189v1.full-text)
<a href="https://doi.org/10.5281/zenodo.16533872"><img src="https://zenodo.org/badge/715874526.svg" alt="DOI"></a>

**IsoLamp is a bash pipeline for the identification of known and novel isoforms from targeted amplicon long-read sequencing data generated with Oxford Nanopore (ONT) sequencing.**

<img src="https://github.com/josiegleeson/IsoLamp/assets/30969357/36e878d7-46c8-4f6a-a2ef-fb39d42e2339" width="960" height="440">

## Contents

- [Installation](#installation)
- [General Command Line Usage](#general-command-line-usage)
- [SIRV Test Dataset](#SIRV-test-dataset)
- [Dependencies](#dependencies)
- [Input Data](#input-data)
- [Parameters](#parameters)
- [Output](#Output)
- [Visualisation with IsoVis](#Visualisation-with-IsoVis)
- [Citation](#Citation)

## Installation
Download the pipeline with Git:
```
git clone https://github.com/ClarkLaboratory/IsoLamp.git
```

Make the script executable and create a conda environment with required dependencies:
```
cd IsoLamp
chmod +x IsoLamp
conda env create -f IsoLamp_env.yml
conda activate IsoLamp
```

## General Command Line Usage
We suggest adding IsoLamp to your PATH so it can be run from any directory:
```
echo 'export PATH=~/path/to/IsoLamp:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

Run IsoLamp (after editing the parameters file):
```
conda activate IsoLamp

# default
IsoLamp -p params.ini
# or
IsoLamp -p params.ini -m all

# for help
IsoLamp -h

# isoform discovery modules only (requires BAMs)
IsoLamp -p params.ini -m isoform_discovery

# alternatively, if not in PATH:
./IsoLamp -p parameters.ini # if running from same directory
~/path/to/IsoLamp/IsoLamp -p parameters.ini # provide full path when running from other directories
```

## SIRV Test Dataset
Test the installation on the provided SIRV data:
```
conda activate IsoLamp
cd IsoLamp # run from within IsoLamp directory for test data
IsoLamp -p sirv_test_data/sirv_params.ini
```
This should produce a folder called 'SIRV5_test' which contains the expected output of the pipeline.

## Dependencies
All dependencies and R libraries are automatically built in the provied conda environment file. Alternatively, the following packages are required.

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
The pipeline is designed to run on barcoded data produced on the Oxford Nanopore platform, typically from amplicon sequencing of a gene of interest. This is typically a directory containing subdirectories titled 'barcode01', 'barcode02' etc. Each barcode directory contains either single or multiple FASTA or FASTQ files of reads. If you have multiple experiments amplifying different genes, run the pipeline once for each gene.

Please see the sirv_test_data/test_fasta directory for an example of the structure required. If you don't have barcoded data, please create a subdirectory with your single sample name and place all read FASTA or FASTQ files inside this directory (nanopore_data/barcode01/reads.fastq, nanopore_data/barcode02/reads.fastq, ...).

**Note:** the 'barcode' name is arbitrary, your subdirectories can be renamed to sample names or any identifiers.

## Parameters
Basic parameters file to run IsoLamp with default settings:
```
OUTPUT_NAME=CLCN3
ENSG_ID=ENSG00000109572
READS="path/to/sample/fa_folders"
GENOME="path/to/genome.fa"
ANNOTATION="path/to/annotation.gtf"
```
Full parameters file (with comments):
```
## Required

## Experiment info ##
OUTPUT_NAME=CLCN3
ENSG_ID=ENSG00000109572

## Path to reads ##
READS="path/to/sample/fa_folders"

## Path to references ##
GENOME="path/to/genome.fa"
ANNOTATION="path/to/annotation.gtf"

## Optional

## Sample grouping ##
# Optional, path to a CSV of sample IDs and group IDs
grouping_data="path/to/group_ids.csv"

## Minimum isoform expression thresholds ##
TPM_minimum=5000
samples_minimum=4 # defualt is one quarter of total number of samples

## Downsampling ##
# downsampling to the same number of reads for every sample, this reduces inter-sample bias
downsampling=TRUE
number_reads_downsample=10000

## Primer site filter ##
# Recommended to increase accuracy and reduce isoforms to only those amplified with primers
primer_site_based_filter=TRUE 
# Primer coordinates for forward and reverse, as bed files
forward_primers="path/to/forward.bed" 
reverse_primers="path/to/reverse.bed" 

## Read accuracy filter ##
extract_high_accuracy_reads=TRUE
minimum_read_accuracy=0.95

## Splice junction accuracy filter ##
extract_high_quality_SJs=FALSE
# minimum junction alignment quality
JAQ=0.9
# window upstream and downstream of splice junction
junction_window=15

## minimap2 ##
# minimap2 -G intron length
max_intron_length=400

## Bambu ##
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
|TPM_minimum|optional|5000|Known and novel isoforms must meet this threshold in at least 'samples_minimum' number of samples.|
|samples_minimum|optional|NA|Known and novel isoforms must meet the 'TPM_minimum' in this many samples. Default value caclulated as: number of samples/4 rounded down to nearest integer. If only 1 sample is provided, the default value will be set to 1.|
|downsampling|optional|TRUE|Whether each sample/barcode should be downsampled to a consistent number of reads.|
|number_reads_downsample|optional|10000|Number of reads to downsample each sample/barcode.|
|primer_site_based_filter|optional|FALSE|Whether to remove isoforms that do not overlap the primers used to perform amplicon sequencing. We **highly** recommend using this option.|
|forward_primers|optional|NULL|BED file of forward primers. Only checked if primer_site_based_filter is TRUE.|
|reverse_primers|optional|NULL|BED file of reverse primers. Only checked if primer_site_based_filter is TRUE.|
|extract_high_accuracy_reads|optional|TRUE|Whether to extract reads with a high accuracy and run Bambu on these reads. We recommend using this option to increase accuracy in calling novel isoforms. All downsampled reads are used for final quantification.|
|minimum_read_accuracy|optional|0.95|The minimum read accuracy, only checked if extract_high_accuracy_reads is TRUE.|
|extract_high_quality_SJs|optional|FALSE|Whether to extract reads with high quality splice junctions, and run Bambu on these reads. This option may increase accuracy in calling novel isoforms. Extracting SJs significantly increases run time and memory usage.|
|JAQ|optional|0.9|The minimum junction alignment quality every junction must meet in order for a read to be included. Only checked if extract_high_quality_SJs is TRUE.|
|junction_window|optional|15|The nt distance upstream and downstream of splice junctions to calculate the JAQ. Only checked if extract_high_quality_SJs is TRUE.|
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
  - quantification of known and novel isoforms (counts, propoprtions, TPMs) as CSVs
  - a PCA plot of samples/barcodes
  - an accuracy plot of all input reads
  - if applicable, a CSV output of t-test results comparing isoform proportions between groups

## Visualisation with IsoVis
The results from the IsoLamp pipeline can be visualised using IsoVis: https://isomix.org/isovis.

The isoform annotation GTF (and optionally isoform counts CSV) can be directly uploaded to IsoVis.

## Citation
If you use IsoLamp, please cite the following paper.

De Paoli-Iseppi et al. (2024) Long-read sequencing reveals the RNA isoform repertoire of neuropsychiatric risk genes in human brain. medRxiv doi: https://doi.org/10.1101/2024.02.22.24303189.









