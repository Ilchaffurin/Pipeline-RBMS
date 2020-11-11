# AUTOMATED PIPELINE FOR ANALYSIS OF GENETIC VARIATIONS OF rRNA

## Table of contents 
1. [Description](#descrp)
2. [Requirements](#req)
3. [Usage](#usage)
4. [Authors](#authors)
5. [Project status](#project)


<a name="descrp"></a> 

## Description

The aim is to analyse the raw data in order to identify the genetic variations from the dataset of lung and breast cancer, in order to better understand the ribosome's heterogeneity and its role in the translational regulation. 

Aim pipeline 

The pipeline uses several software like : 
* FASTQC for quality control 
* TRIMMOMATIC for trimming and cleaning reads
* STAR for indexing and mapping 
* BOWTIE2 for indexing and mapping 
* BCFTOOLS Identification of variants and creation of vcf files 
* SAMTOOLS Convert a SAM file into a BAM file 
* MULTIQC Aggregate results from all analyses into a single report

Pipeline steps  

The pipeline is organised as follows  (imagen)  

<a name="req"></a> 

## Requirements 

##### Install miniconda3

[Miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers)

##### Create the environment from EnvPipeline.yml file : 
``` conda env create -f EnvPipeline.yml ```

<a name="usage"></a> 

## Usage 

##### To run the pipeline :

```
nextflow run pipeline.nf --input_dir 'path/to/fasta/*.fastq' --genome_ref 'path/to/genome/file.fasta'
```

##### To run the pipeline help message : 

```
nextflow run pipeline.nf --help 
```

##### Help message : 
```
Required arguments:
    --input_dir      Directory for fastq files
    --genome_ref     Full path to directory containing reference genome fasta file

QC option:
    --skipFastqc     Skip reads quality control step (default: activated).
    --skipMultiqc    Skip merging tools reports suitable with multiqc (default: activated)

Trimming option:
    --skipTrmming    Skip trimming step (default: activated).

Mapping option:
    --onlySTAR       Only using STAR mapper (default: STAR and Bowtie2).
    --onlyBowtie2    Only using Bowtie2 mapper (default: STAR and Bowtie2).

Save option:
    --outdir         Specify where to save the output from the nextflow run (default: "./results/")

Threading option   
    --threads        Specify the number of threads to be used (default: 1).

help message:
    --help           Print help message
```



<a name="authors"></a> 

## Authors 

The project was developped by : 

* Chloé AUJOULAT
* Ilan CHAFFURIN
* Morgane DES LIGNERIS
* Isis LORENZO-COLINA


<a name="project"></a> 

## Project status 

