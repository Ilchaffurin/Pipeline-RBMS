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
* `FASTQC` for quality control 
* `TRIMMOMATIC` for trimming and cleaning reads
* `STAR` for indexing and mapping 
* `BOWTIE2` for indexing and mapping 
* `BCFTOOLS` to identify variants and creation of cvs files 
* `SAMTOOLS` to convert and sort a SAM file into a BAM file 
* `MULTIQC` Aggregate results from all analyses into a single report

#### The pipeline is organized as follows:  

![alt text](/img/schema-pipeline.png)

<a name="req"></a> 

## Requirements 

#### Install miniconda3: 

[Miniconda](https://docs.conda.io/en/latest/miniconda.html#linux-installers)

#### Create the environment from EnvPipeline.yml file: 
``` conda env create -f EnvPipeline.yml ```

#### Activate the environement EnvPipeline:
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
    --dp             Read length option (default: 50)

Miniconda3 path option:
    --miniconda3     Full path to directory containing miniconda3 (exemple : '/data/home/.../miniconda3')

QC option:
    --skipFastqc     Skip reads quality control step (default: activated).
    --skipMultiqc    Skip merging tools reports suitable with multiqc (default: activated)

Trimming option:
    --skipTrimming   Skip trimming step (default: activated).

Mapping option:
    --onlySTAR       Only using STAR mapper (default: STAR and Bowtie2).
    --onlyBowtie2    Only using Bowtie2 mapper (default: STAR and Bowtie2).
    --skipMapping    Skip mapping step (default: activated).

Save option:
    --outdir         Specify where to save the output from the nextflow run (default: "./results/")

Threading option   
    --threads        Specify the number of threads to be used (default: 1).

help message:
    --help           Print help message
```

#### To run the pipeline:

```nextflow
nextflow run pipeline.nf --threads 10 --input_dir 'path/to/*.fastq' --genome_ref 'path/to/file.fasta' --miniconda3 'path/to/miniconda3'
```

#### Execution example: 

```nextflow
  nextflow run pipeline.nf --threads 15 --input_dir '/data/home/mdes-ligneris/M1/S1/Projet1/data/B*.fastq' --genome_ref '/data/home/mdes-ligneris/M1/S1/Projet1/data/Human.fa' --miniconda3 '/data/home/mdes-ligneris/miniconda3' --bp 75
 ```

#### After execution the script will output all the information process : 

```shell
N E X T F L O W  ~  version 20.10.0
Launching `pipeline.nf` [awesome_kilby] - revision: 12da2a83da
------------------------------------------------------------------------------------------------------------------------
Fastq file(s) from Path  : /data/home/mdes-ligneris/M1/S1/Projet1/data/B*.fastq
Genome fasta file        : /data/home/mdes-ligneris/M1/S1/Projet1/data/Human.fa
Output                   : results
Number of threads        : 15
Fragment length          : 75
Reads QC                 : Yes
Merging Reports          : Yes
Trimming                 : Yes
Mapper                   : STAR and Bowtie2
------------------------------------------------------------------------------------------------------------------------
executor >  local (45)
[8b/acc838] process > Trimmomatic (B3199_10_S8_R1_001) [100%] 6 of 6 ✔
[a5/b13d01] process > Fastqc (B3199_6_S7_R1_001)       [100%] 12 of 12 ✔
[d6/6b19d5] process > Index_STAR (Human)               [100%] 1 of 1 ✔
[d9/627d77] process > Index_BOWTIE (Human)             [100%] 1 of 1 ✔
[80/7f03e2] process > Star (B3199_6_S7_R1_001)         [100%] 6 of 6 ✔
[f4/f7bb8b] process > Bowtie2 (B3446_6_S9_R1_001)      [100%] 6 of 6 ✔
[72/bfacac] process > Bcftools (B3446_6_S9_R1_001)     [100%] 12 of 12 ✔
[69/cb74f1] process > MultiQC                          [100%] 1 of 1 ✔
Completed at: 15-Dec-2020 19:44:22
Duration    : 49m 29s
CPU hours   : 11.2
Succeeded   : 45
```

#### Where are the outputs files ? 

By default, all the files are saved in the `./results` directory. In the directory you will find five different directories named as each process:  

* FastQC for the HTML quatity control files before and after trimming and cleaning 
* Trimming for the fastq file(s) after trimming and cleaning  
* mapping for the index and align directory for STAR and BOWTIE2 process 
* multiQC for the HTML file for the final quality rapport 
* bcftools with the vcf and cvs files for the variant calling process 

#### Convert VCF files into CVS files 

If you want to change the files separators, you can modify the `;` (semicolon) in the string and insert the separator of interest. 

The command concerned can be found in the `variant calling` part with the name `process Bcftools` 

* In the following example the separators `;` (semicolon) were replaced by `,` (comma)

```
bcftools query -f '%CHROM,%POS,%ID,%REF,%ALT,%QUAL,%FILTER\n' ${file_id}${data_type}calls.vfc -o ${file_id}${data_type}_calls.csv -H
```

#### 

#### Number of mismatches 

By default, the number of mismatch autorized by the mappers is one. To change this parameter you need to replace the '1' by the number of mismatched wanted in ligne 130 for STAR parameters :
```
129  if (params.bp){
130    lmax = 1/(params.bp.toInteger())
131  }
```

and the '1' by the number of mismatched wanted behind '-N' in ligne 444 for bowtie2 :
```
443  bowtie2 -p ${cpus} \
444            -D 20 -R 2 -N 1 -L 20 -i S,1,0.50 \
445            -x ${index_id} \
446            -U ${reads} \
447            -S ${file_id}_bowtie2_tmp.sam 2> ${file_id}_bowtie2.stats.log
```

#### 

<a name="authors"></a> 

## Authors 

The project was developped by : 

* Chloé AUJOULAT
* Ilan CHAFFURIN
* Morgane DES LIGNERIS
* Isis LORENZO-COLINA


<a name="project"></a> 



