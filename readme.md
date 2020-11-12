# AUTOMATED PIPELINE FOR ANALYSIS OF GENETIC VARIATIONS OF rRNA

## Table of contents 
1. [Description](#descrp)
2. [Requirements](#req)
3. [Usage](#usage)
4. [Authors](#authors)
5. [Project status](#project)


<a name="descrp"></a> 

## Description

The aim of this project is to analyse the fastq raw data in order to identify the genetic variations from the dataset of lung and breast cancer, with the pipeline we will be able to better understand the ribosome's heterogeneity and its role in the translational regulation. 

For the analyse we establish a scalable and reproducible pipeline in Nextflow workflow. 

The pipeline uses several software like : 
* FASTQC for quality control 
* TRIMMOMATIC for trimming and cleaning reads
* STAR for indexing and mapping 
* BOWTIE2 for indexing and mapping 
* BCFTOOLS to identify variants and creation of cvs files 
* SAMTOOLS to convert a SAM file into a BAM file 
* MULTIQC Aggregate results from all analyses into a single report

##### The pipeline is organized as follows:  

![alt text](/img/pipeline.png)

<a name="req"></a> 

## Requirements 

##### Install miniconda3

[Miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers)

##### Create the environment from EnvPipeline.yml file: 
``` conda env create -f EnvPipeline.yml ```

##### Activate the environement EnvPipeline 
``` conda activate EnvPipeline ```

<a name="usage"></a> 

## Usage 

#### To run the pipeline help message: 

```nextflow
nextflow run pipeline.nf --help 
```

#### Help message: 
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

#### To run the pipeline:

```nextflow
nextflow run pipeline.nf --threads 10 --input_dir 'path/to/*.fastq' --genome_ref 'path/to/file.fasta'
```

#### Execution example: 

```nextflow
 nextflow run pipeline.nf --input_dir '/data/home/mdes-ligneris/M1/S1/Projet1/data/B2998_10_S6_R1_001.fastq' --genome_ref /data/home/mdes-ligneris/M1/S1/Projet1/data/Human.fa --outdir resultsCBoth
 ```

#### After execution the script will output all the process : 

```
N E X T F L O W  ~  version 20.10.0
Launching `pipeline.nf` [gloomy_linnaeus] - revision: d53c8822c5
------------------------------------------------------------------------------------------------------------------------
Fastq file(s) from Path  : /data/home/mdes-ligneris/M1/S1/Projet1/data/B2998_10_S6_R1_001.fastq
Genome fasta file        : /data/home/mdes-ligneris/M1/S1/Projet1/data/Human.fa
Output                   : results
Number of threads        : 1
Reads QC                 : Yes
Merging Reports          : Yes
Trimming                 : Yes
Mapper                   : STAR ans Bowtie2
------------------------------------------------------------------------------------------------------------------------
executor >  local (7)
[6d/070e25] process > Trimmomatic (B2998_10_S6_R1_001) [100%] 1 of 1 ✔
[0e/e6067a] process > Fastqc (B2998_10_S6_R1_001)      [100%] 2 of 2 ✔
executor >  local (7)
[6d/070e25] process > Trimmomatic (B2998_10_S6_R1_001) [100%] 1 of 1 ✔
[0e/e6067a] process > Fastqc (B2998_10_S6_R1_001)      [100%] 2 of 2 ✔
[a8/d9bdab] process > Index_STAR (Human)               [100%] 1 of 1 ✔
[f2/f96b04] process > Index_BOWTIE (Human)             [100%] 1 of 1 ✔
[-        ] process > Star (B2998_10_S6_R1_001)        -
[fe/d28b0e] process > Bowtie2 (B2998_10_S6_R1_001)     [100%] 1 of 1 ✔
[-        ] process > Bcftools                         -
[-        ] process > MultiQC                          -
```

#### Where are the outputs files ? 

By default, all the files are saved in the `./results` directory. In the directory you will find five different file with the name of each process ()

#### Convert VCF files into CVS files 

If you want to change the files separators, you can modify the `;` in the string and insert the separator of interest 

Example: 

* In the script the separators are indicated by `;` (semicolon)

```
bcftools query -f '%CHROM;%POS;%ID;%REF;%ALT;%QUAL;%FILTER\n' ${file_id}${data_type}calls.vfc -o ${file_id}${data_type}_calls.csv -H
```

* In the following example the separators are indicated by `,` (comma)

```
bcftools query -f '%CHROM,%POS,%ID,%REF,%ALT,%QUAL,%FILTER\n' ${file_id}${data_type}calls.vfc -o ${file_id}${data_type}_calls.csv -H
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


