#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      nextflow run pipeline.nf --fasta <fasta files> 

    Required arguments:
      --fastq           Directory pattern for fastq files

    Save option:
      --outdir          Specify where to save the output from the nextflow run (default: "./results/")

    help message:
      --help            Print help message
    """
      .stripIndent()
  }

///////////////////////////////////////////////////////////////////////////////
/* --                   DEFAULT PARAMETER VALUES                          -- */
///////////////////////////////////////////////////////////////////////////////

params.help = false
params.fastq = false
params.outdir = 'results'

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

///////////////////////////////////////////////////////////////////////////////
/* --                       HEADER LOG INFO                               -- */
///////////////////////////////////////////////////////////////////////////////

//log.info "fastq files    : ${params.fastq}"

def summary = [:]
summary['Fastq files']              = params.fastq ? params.fastq : 'Not supplied'
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")

///////////////////////////////////////////////////////////////////////////////
/* --                           VALIDATE INPUTS                           -- */
///////////////////////////////////////////////////////////////////////////////

//params.fastq = "/data/home/mdes-ligneris/M1/S1/Projet1/data/B2998_10_S6_R1_001.fastq"

if (params.fastq){
    Channel
        .fromFilePairs( params.fastq,size:1 )
        .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq}" }
        .into { fastqc_files;fastq_files_for_report }
}
else { exit 1,
    log.warn "=================================================================\nWARNING! No fastq files precised.\nUse '--fastq' \nOr '--help' for more informations"
}

///////////////////////////////////////////////////////////////////////////////
/* --                       READS QUALITY CONTROLE                        -- */
///////////////////////////////////////////////////////////////////////////////

process Fastqc {
    label "fastqc"
        tag "$file_id"
        publishDir "${params.outdir}/fastq/QC/", mode: 'copy'
          
        input:
        set file_id, file(reads) from fastq_files_for_report

        output:
        file "*.{zip,html}" into fastqc_report

        script:
        """
        fastqc ${reads} --format fastq --outdir ./
        """
}

 ///////////////////////////////////////////////////////////////////////////////
 /* --                       TRIMMING / CLEANING READS                     -- */
 ///////////////////////////////////////////////////////////////////////////////

