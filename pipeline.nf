#!/usr/bin/env nextflow

params.fastq = "/data/home/mdes-ligneris/M1/S1/Projet1/data/B2998_10_S6_R1_001.fastq"

params.outdir = 'results'

log.info "fastq files : ${params.fastq}"

Channel
  .fromFilePairs( params.fastq,size:1 )
  .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
  .set { fastq_files }

process Fastqc {
    label "fastqc"
        tag "$file_id"
        publishDir "${params.outdir}/fastq/QC/", mode: 'copy'
          
        input:
        set file_id, file(reads) from fastq_files

        output:
        file "*.{zip,html}" into fastqc_report

        script:
        """
        fastqc ${reads} --format fastq --outdir ./

        """
}
