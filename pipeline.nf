#!/usr/bin/env nextflow

params.fastq = "mdes-ligneris@pedago-ngs:~/M1/S1/Projet1/data/B2998_10_S6_R1_001.fastq"

params.outdir = 'results'

log.info "fastq files : ${params.fastq}"

Channel
  .fromPath( params.fastq )
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
        fastqc --threads ${reads} --format fastq --outdir ./
        """
}