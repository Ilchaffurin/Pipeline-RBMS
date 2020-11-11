# AUTOMATED PIPELINE FOR ANALYSIS OF GENETIC VARIATIONS OF rRNA

## Table of contents 
1. [Description](#descrp)
2. [Requirements](#req)
3. [Usage](#usage)
4. [Authors and ackowledment](#authors)
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
* MULTIQC 

Pipeline steps  

The pipeline is organised as follows  (imagen)  

<a name="req"></a> 

## Requirements 

##### Install miniconda3

https://docs.conda.io/en/latest/miniconda.html#linux-installers - automatic!


##### Create the environment from EnvPipeline.yml file : 
``` conda env create -f EnvPipeline.yml ```

<a name="usage"></a> 

## Usage 

<a name="authors"></a> 

## Authors and ackowledment 

The project was developped by : 

* Chloé AUJOULAT
* Ilan CHAFFURIN
* Morgane DES LIGNERIS
* Isis LORENZO-COLINA


<a name="project"></a> 

## Project status 

