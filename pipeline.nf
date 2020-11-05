#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      nextflow run pipeline.nf --fasta <fasta files> 

    Required arguments:
      --input_dir      Directory pattern for fastq files

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
params.input_dir = false
params.outdir = 'results'

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

///////////////////////////////////////////////////////////////////////////////
/* --                          HEADER LOG INFO                            -- */
///////////////////////////////////////////////////////////////////////////////

log.info "fastq files    : ${params.input_dir}"

///////////////////////////////////////////////////////////////////////////////
/* --                          VALIDATE INPUTS                            -- */
///////////////////////////////////////////////////////////////////////////////

if (params.input_dir){
    Channel
        .fromFilePairs( params.input_dir , size:1 )
        .ifEmpty { error "Cannot find any fastq files matching: ${params.input_dir}" }
        .into { fastqc_files_2trim ; fastq_files_2QC }
}
else { 
    log.info "No fastq files precised.\nUse '--fastq' \nOr '--help' for more informations"
    exit 1
}

///////////////////////////////////////////////////////////////////////////////
/* --                     FIRST READS QUALITY CONTROLE                    -- */
///////////////////////////////////////////////////////////////////////////////

process Fastqc {
    label "fastqc"
        tag "$file_id"
        publishDir "${params.outdir}/fastq/QC/", mode: 'copy'
          
        input:
        set file_id, file(reads) from fastq_files_2QC

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
        trimmomatic SE -phred33 ${reads} ${file_id}_trim ILLUMINACLIP:trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:7 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:30 MINLEN:8
        """
}

