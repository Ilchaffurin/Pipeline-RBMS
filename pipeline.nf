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

    QC Option:
      --skipFastqc                  Skip reads quality control step (default: activated).
      --skipMultiqc                 Skip merging tools reports suitable with multiqc (default: activated)


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
params.skipFastqc = false
params.skipMultiqc = false
params.outdir = 'results'

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

///////////////////////////////////////////////////////////////////////////////
/* --                          HEADER LOG INFO                            -- */
///////////////////////////////////////////////////////////////////////////////

log. info "-------------------------------------------------------------------------------"
log.info "path to fastq files  : ${params.input_dir}"
log.info "genome fasta file    : ${params.index}"
log.info "output               : ${params.outdir}"
if (params.skipFastqc){
  log.info "Reads QC             : Skipped"
}
else {
  log.info "Reads QC             : Yes"
}
if (params.skipMultiqc){
  log.info "Merging Reports      : Skipped"
}
else {
  log.info "Merging Reports      : Yes"
}
log. info "-------------------------------------------------------------------------------"

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
        .into { fasta_2indexing_STAR ; fasta_2indexing_BOWTIE; fasta_2variantCalling }
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
        set file_id, "*.fastq" into fastq_trim_files, fastq_trim_files2, fastq_trim_files_2QC, trimming_report

        script:
        """
        trimmomatic SE -phred33 ${reads} ${file_id}_trim.fastq \
        ILLUMINACLIP:~/miniconda2/envs/EnvPipeline/share/trimmomatic-0.39-1/adapters/TruSeq3-SE.fa:2:30:7 \
        LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:30 MINLEN:8
        """
}

///////////////////////////////////////////////////////////////////////////////
/* --                         READS QUALITY CONTROLE                      -- */
///////////////////////////////////////////////////////////////////////////////

fastq_files_2QC
        .concat(fastq_trim_files_2QC)
        .set { fastq_files }

if (params.skipFastqc) {
         Channel
               .empty()
               .set { fastqc_report }
 }
 else{
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
}

///////////////////////////////////////////////////////////////////////////////
/* --                               INDEX                                 -- */
///////////////////////////////////////////////////////////////////////////////

process Index_STAR {
    label "index_star"
        tag "$file_id"
        publishDir "${params.outdir}/mapping/STAR/index/", mode: 'copy'
          
        input:
        set file_id, file(fasta) from fasta_2indexing_STAR

        output:
        file "*" into index_files_STAR

        script:
        """
        STAR --runThreadN 20 \
        --runMode genomeGenerate \
        --genomeDir ./index/ \
        --genomeFastaFiles ${fasta} 
        """
}

process Index_BOWTIE {
    label "index_bowtie2"
        tag "$file_id"
        publishDir "${params.outdir}/mapping/BOWTIE2/index/", mode: 'copy'
          
        input:
        set file_id, file(fasta) from fasta_2indexing_BOWTIE

        output:
        file "*.index*" into index_files_BOWTIE
	      file "*_report.txt" into indexing_report

        script:
        """
        bowtie2-build ${fasta} ${file_id}.index &> ${fasta.baseName}_bowtie2_report.txt
        """
}


///////////////////////////////////////////////////////////////////////////////
/* --                              MAPPING                                -- */
///////////////////////////////////////////////////////////////////////////////

process Star {
    label "star"
        tag "$file_id"
        publishDir "${params.outdir}/mapping/STAR/$file_id", mode: 'copy'
          
        input:
        set file_id, file(reads) from fastq_trim_files
	      file index from index_files_STAR.collect()

        output:
        set file_id, val(data_type), "*out.bam" into star_bam_files
        file "*" into star_mapping_report

        script:
        data_type="_star_"
        """
        STAR --runThreadN 20 \
        --runMode alignReads \
        --genomeDir ./index/ \
        --readFilesIn ${reads} \
        --outFileNamePrefix ./${file_id} \
        --outSAMtype BAM SortedByCoordinate
        """
}

process Bowtie2 {
    label "bowtie2"
        tag "$file_id"
        publishDir "${params.outdir}/mapping/BOWTIE2/$file_id", mode: 'copy'
          
        input:
        set file_id, file(reads) from fastq_trim_files2
	      file index from index_files_BOWTIE.collect()

        output:
        set file_id, val(data_type), "*_sorted.bam" into bowtie2_bam_files
        file "*" into bowtie2_mapping_report

        script:
        data_type="_bowtie2_"
        index_id = index[0]
        for (index_file in index) {
          if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
            index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
          }
        }
        """
        bowtie2 -p 20 \
        --very-sensitive \
        -x ${index_id} \
        -U ${reads} \
        -S ${file_id}_bowtie2.sam 
        samtools view -bS ${file_id}_bowtie2.sam > ${file_id}_bowtie2.bam
        samtools sort -@ 20 -O BAM -o ${file_id}_bowtie2_sorted.bam ${file_id}_bowtie2.bam
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
        set file_id, data_type, file(reads) from star_bam_files
        set file_id2, data_type2, file(reads2) from bowtie2_bam_files
	      set file_id3, file(fasta) from fasta_2variantCalling
        
        output:
        set file_id, "*" into variant_calling_file

        script:
        """
        bcftools mpileup -f $fasta $reads | bcftools call -mv -Ob -o ${file_id}${data_type}calls.vcf
        bcftools mpileup -f $fasta $reads2 | bcftools call -mv -Ob -o ${file_id2}${data_type2}calls.vcf

	      bcftools view -i '%QUAL>=20' ${file_id}${data_type}calls.vcf > ${file_id}${data_type}calls_view.vcf
        bcftools view -i '%QUAL>=20' ${file_id2}${data_type2}calls.vcf > ${file_id2}${data_type2}calls_view.vcf

        """     
}

///////////////////////////////////////////////////////////////////////////////
/* --                      MERGE ALL STEPS REPORTS                        -- */
///////////////////////////////////////////////////////////////////////////////

process MultiQC {
    label "multiQC"
      publishDir "${params.outdir}/multiQC", mode: 'copy'

      input:
      file report_fastqc from fastqc_report.collect().ifEmpty([])
      file report_trim from trimming_report.collect().ifEmpty([])
      file report_star from star_mapping_report.collect().ifEmpty([])
      file report_mapping from bowtie2_mapping_report.collect()

      output:
      file "*multiqc_*" into multiqc_report

      when:
      !params.skipMultiqc

      script:
      """
      multiqc -f . \\
      -m fastqc -m trimmomatic -m star -m bowtie2
      """
}
