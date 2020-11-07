#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      nextflow run pipeline.nf --fasta <fasta files> 

    Required arguments:
      --input_dir      Directory for fastq files
    
    Reference genome
      --index          Full path to directory containing genome fasta file

    Save option:
      --outdir         Specify where to save the output from the nextflow run (default: "./results/")

    help message:
      --help           Print help message
    """
      .stripIndent()
  }

///////////////////////////////////////////////////////////////////////////////
/* --                   DEFAULT PARAMETER VALUES                          -- */
///////////////////////////////////////////////////////////////////////////////

params.help = false
params.input_dir = false
params.index = false
params.outdir = 'results'

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

///////////////////////////////////////////////////////////////////////////////
/* --                          HEADER LOG INFO                            -- */
///////////////////////////////////////////////////////////////////////////////

log.info "fastq files    :       ${params.input_dir}"
log.info "genome fasta file    : ${params.index}"

///////////////////////////////////////////////////////////////////////////////
/* --                          VALIDATE INPUTS                            -- */
///////////////////////////////////////////////////////////////////////////////

if (params.input_dir){
    Channel
        .fromFilePairs( params.input_dir , size:-1 )
        .ifEmpty { error "Cannot find any fastq files matching: ${params.input_dir}" }
        .into { fastqc_files_2trim ; fastq_files_2QC }
}
else { 
    log.info "No fastq files precised.\nUse '--input_dir' \nOr '--help' for more informations"
    exit 1
}

if (params.index){
    Channel
        .fromFilePairs( params.index , size:1 )
        .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
        .into { fasta_2indexing ; fasta_2variantCalling }
}
else { 
    log.info "No index file precised.\nUse '--index' \nOr '--help' for more informations"
    exit 1
}

///////////////////////////////////////////////////////////////////////////////
/* --                       TRIMMING / CLEANING READS                     -- */
///////////////////////////////////////////////////////////////////////////////

process Trimmomatic {
    label "trimmomatic"
        tag "$file_id"
        publishDir "${params.outdir}/fastq/trim/", mode: 'copy'
          
        input:
        set file_id, file(reads) from fastqc_files_2trim

        output:
        set file_id, "*.fastq" into fastq_trim_files, fastq_trim_files_2QC

        script:
        """
        trimmomatic SE -phred33 ${reads} ${file_id}_trim.fastq ILLUMINACLIP:~/miniconda2/envs/EnvPipeline/share/trimmomatic-0.39-1/adapters/TruSeq3-SE.fa:2:30:7 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:30 MINLEN:8
        """
}

///////////////////////////////////////////////////////////////////////////////
/* --                         READS QUALITY CONTROLE                      -- */
///////////////////////////////////////////////////////////////////////////////

process Fastqc {
    label "fastqc"
        tag "$file_id"
        publishDir "${params.outdir}/fastq/QC/", mode: 'copy'
          
        input:
        set file_id, file(reads) from fastq_files_2QC
        set file_id, file(reads) from fastq_trim_files_2QC

        output:
        file "*.{zip,html}" into fastqc_report

        script:
        """
        fastqc ${reads} --format fastq --outdir ./
        """
}

///////////////////////////////////////////////////////////////////////////////
/* --                               INDEX                                 -- */
///////////////////////////////////////////////////////////////////////////////

process Index_STAR {
    label "index_star"
        tag "$file_id"
        publishDir "${params.outdir}/mapping/STAR/index/", mode: 'copy'
          
        input:
        set file_id, file(reads) from fasta_2indexing

        output:
        file "*" into index_files

        script:
        """
        STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ./index/ --genomeFastaFiles ${reads} 

        """
}

///////////////////////////////////////////////////////////////////////////////
/* --                              MAPPING                                -- */
///////////////////////////////////////////////////////////////////////////////

process STAR {
    label "star"
        tag "$file_id"
        publishDir "${params.outdir}/mapping/STAR/$file_id", mode: 'copy'
          
        input:
        set file_id, file(reads) from fastq_trim_files
	file index from index_files.collect()

        output:
        set file_id, "*out.bam" into bam_file
        file "*" into star_report

        script:
        """
        STAR --runThreadN 10 \
        --runMode alignReads \
        --genomeDir ./index/ \
        --readFilesIn ${reads} \
        --outFileNamePrefix ./${file_id} \
        --outSAMtype BAM SortedByCoordinate
        """
}

///////////////////////////////////////////////////////////////////////////////
/* --                         VARIANT CALLING                             -- */
///////////////////////////////////////////////////////////////////////////////

process Bcftools {
    label "bcftools"
        tag "$file_id"
        publishDir "${params.outdir}/bcftools/$file_id", mode: 'copy'
          
        input:
        set file_id, file(reads) from bam_file
	set file_id2, file(fasta) from fasta_2variantCalling

        output:
        set file_id, "*" into variant_calling_file

        script:
        """
        bcftools mpileup -f $fasta $reads | bcftools call -mv -Ob -o calls.vcf
	bcftools view -i '%QUAL>=20' calls.vcf > calls_view.vcf
 
        """
}
