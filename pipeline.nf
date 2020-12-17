#!/usr/bin/env nextflow

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        AUTOMATED PIPELINE FOR                       -- */
/* --               ANALYSIS OF GENETIC VARIATIONS OF rRNA                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/* 
* Setting the User Help Message
*/
def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      nextflow run pipeline.nf --threads 10 --input_dir 'path/to/*.fastq' --genome_ref 'path/to/file.fasta' --miniconda3 'path/to/miniconda3'

    Required arguments:
      --input_dir      Directory for fastq files
      --genome_ref     Full path to directory containing reference genome fasta file
      --dp             Read length option (default: 50)

    Miniconda3 path option:
      --miniconda3     Full path to directory containing miniconda3 (exemple : '/data/home/.../miniconda3')

    QC option:
      --skipFastqc     Skip reads quality control step (default: not activated).
      --skipMultiqc    Skip merging tools reports suitable with multiqc (default: not activated)

    Trimming option:
      --skipTrimming   Skip trimming step (default: not activated).

    Mapping option:
      --onlySTAR       Only using STAR mapper (default: STAR and Bowtie2).
      --onlyBowtie2    Only using Bowtie2 mapper (default: STAR and Bowtie2).
      --skipMapping    Skip mapping step (default: not activated).

    Save option:
      --outdir         Specify where to save the output from the nextflow run (default: "./results/")

    Threading option   
      --threads        Specify the number of threads to be used (default: 1).

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
params.genome_ref = false
params.skipFastqc = false
params.skipMultiqc = false
params.skipTrimming = false
params.skipMapping = false
params.onlySTAR = false
params.onlyBowtie2 = false
params.outdir = 'results'
params.miniconda3= '~/miniconda3'
params.threads = 1
params.bp = '50'
max50 = 0.04
max75 = 0.027

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

///////////////////////////////////////////////////////////////////////////////
/* --                          VALIDATE INPUTS                            -- */
///////////////////////////////////////////////////////////////////////////////

/* 
* Initialization of all fastq files returned by the path into tuples in which 
* the first element is the key of the file and the second element is the file. 
*/
if (params.input_dir){
    Channel
        .fromFilePairs(params.input_dir, size:1, checkIfExists: true)
        .ifEmpty { error "Cannot find any fastq files matching: ${params.input_dir}" }
        .into { fastqc_files_2trim ; fastq_files_2QC }    
}
else { 
    log.info "No fastq files precised.\nUse '--input_dir' \nOr '--help' for more informations"
    exit 1
}

/* 
* Initialization of reference genome fasta file returned by the path into tuples 
* in which the first element is the ID of the file and the second element is the 
* file. The ID is created from the name of the file without the extension.
*/
if (params.genome_ref){
    Channel
        .fromFilePairs( params.genome_ref , size:1, checkIfExists: true )
        .ifEmpty { error "Cannot find any genome_ref files matching: ${params.genome_ref}" }
        .into { fasta_2indexing_STAR ; fasta_2indexing_BOWTIE; fasta_2variantCalling }
}
else { 
    log.info "No genome_ref file precised.\nUse '--genome_ref' \nOr '--help' for more informations"
    exit 1
}

/* 
* Initialization of the number of threads to be used. 
*/
if (params.threads){
    cpus = "${params.threads}"
}

/* 
* Initialization of the max number of mismatches relative to the read length to be used.
  Needs to be set for the option --outFilterMismatchNoverLmax in Star.
        if 75bp and 1 mismatch max : 1/75 = 0.013 
        if 50bp and 1 musmatch max : 1/50 = 0.02
  The number 1 bellow need to be changed according to the number of mismatch wanted. 
*/
if (params.bp){
  lmax = 1/(params.bp.toInteger())
}

///////////////////////////////////////////////////////////////////////////////
/* --                          HEADER LOG INFO                            -- */
///////////////////////////////////////////////////////////////////////////////

/* 
* Writes all information according to the options at the beginning of the run.
*/
log.info "".padRight(120,'-')
log.info "Fastq file(s) from Path  : ${params.input_dir}"
log.info "Genome fasta file        : ${params.genome_ref}"
log.info "Output                   : ${params.outdir}"
log.info "Number of threads        : ${params.threads}"
log.info "Fragment length          : ${params.bp} "
if (params.skipFastqc){ 
  log.info "Reads QC                 : Skipped"
} 
else {
  log.info "Reads QC                 : Yes"
}
if (params.skipMultiqc){ 
  log.info "Merging Reports          : Skipped"
}
else {
  log.info "Merging Reports          : Yes"
}
if (params.skipTrimming){
  log.info "Trimming                 : Skipped"
}
else {
  log.info "Trimming                 : Yes"
}
if (!params.onlySTAR && !params.onlyBowtie2 && !params.skipMapping){
  log.info "Mapper                   : STAR and Bowtie2"
}
if (params.onlySTAR){
  log.info "Mapper                   : STAR"
}
if (params.onlyBowtie2){
  log.info "Mapper                   : Bowtie2"
}
if (params.skipMapping){
  log.info "Mapper                   : Skipped"
}
log.info "".padRight(120,'-')


///////////////////////////////////////////////////////////////////////////////
/* --                       TRIMMING / CLEANING READS                     -- */
///////////////////////////////////////////////////////////////////////////////

/* 
* Initialization of default empty values if the trimming process is skipped.
*/
if (params.skipTrimming) {
    fastqc_files_2trim.into{ fastq_trim_files ; fastq_trim_files2 }
    Channel
        .empty()
        .into { fastq_trim_files_2QC ; trimming_report }
}

/* 
* Cleaning and trimming of the reads from the fastq file(s). The results are 
* published inside a directory name after the ID of the file.
*/
else {
    process Trimmomatic {
        label "trimmomatic"
        tag "$file_id"
        publishDir "${params.outdir}/Trimming/$file_id", mode: 'copy'
              
        input:
        set file_id, file(reads) from fastqc_files_2trim

        output:
        set file_id, "*.fastq" into fastq_trim_files, fastq_trim_files2, fastq_trim_files_2QC
        file "*trim.fastq.stats.log" into trimming_report

        script:
        """
        trimmomatic SE -threads ${cpus} -phred33 ${reads} ${file_id}_trim.fastq \
        ILLUMINACLIP:${params.miniconda3}/envs/EnvPipeline/share/trimmomatic-0.39-1/adapters/TruSeq3-SE.fa:2:30:7 \
        LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:30 MINLEN:8 2> ${file_id}_trim.fastq.stats.log
        """
    }
}

///////////////////////////////////////////////////////////////////////////////
/* --                         READS QUALITY CONTROLE                      -- */
///////////////////////////////////////////////////////////////////////////////

/* 
* Initialization of default empty values if the FastQC process is skipped.
*/
if (params.skipFastqc) {
         Channel
               .empty()
               .set { fastqc_report }
}
/* 
* Concat the fastq file from before the trimming process, and after it into a 
* single Channel
*/
else {
    fastq_files_2QC
        .concat(fastq_trim_files_2QC)
        .set { fastq_files }
/* 
* Quality evaluation of the fastq files with FastQC. The results are published 
* inside a directory name after the ID of the file.
*/
    process Fastqc {
        label "fastqc"
        tag "$file_id"
        publishDir "${params.outdir}/FastQC/$file_id", mode: 'copy'
          
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

/* 
* Initialization of default empty values if the only the mapper Bowtie2 is used.
*/
if (params.onlyBowtie2 || params.skipMapping){
    Channel
        .empty()
        .set { fasta_2indexing_STAR }
}
/* 
* Creating a reference genome index from the fasta file for STAR mapping.
*/
else {
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
        STAR --runThreadN ${cpus} \
        --runMode genomeGenerate \
        --genomeDir ./index/ \
        --genomeFastaFiles ${fasta} 
        """
  }  
}

/* 
* Initialization of default empty values if the only the mapper STAR is used.
*/
if (params.onlySTAR || params.skipMapping){
    Channel
        .empty()
        .set { fasta_2indexing_BOWTIE }
}
/* 
* Creating a reference genome index from the fasta file for STAR mapping.
*/
else {
    process Index_BOWTIE {
        label "index_bowtie2"
            tag "$file_id"
            publishDir "${params.outdir}/mapping/BOWTIE2/index/", mode: 'copy'
            
            input:
            set file_id, file(fasta) from fasta_2indexing_BOWTIE

            output:
            file "*.index*" into index_files_BOWTIE
            file "*stats.log" into indexing_report

            script:
            """
            bowtie2-build ${fasta} ${file_id}.index 2> ${fasta.baseName}_bowtie2_stats.log
            """
  }
}

///////////////////////////////////////////////////////////////////////////////
/* --                              MAPPING                                -- */
///////////////////////////////////////////////////////////////////////////////

/* 
* Initialization of default empty values if the only the mapper Bowtie2 is used.
*/
if (params.onlyBowtie2 || params.skipMapping){
  Channel
      .empty()
      .into { star_bam_files ; star_mapping_report}
}
/* 
* Alignment of reads using STAR mapper.
*/
else {
    process Star {
        label "star"
        tag "$file_id"
        publishDir "${params.outdir}/mapping/STAR/$file_id", mode: 'copy'
            
        input:
        set file_id, file(reads) from fastq_trim_files
        file index from index_files_STAR.collect()

        output:
        set file_id, val(data_type), "*out.bam" into star_bam_files
        file "*Log.final.out" into star_mapping_report 

        script:
        data_type="_star_"
        
        /*
        Star takes a STAR index and a set of sequencing read files as input (files in fastQ format) and generates a set of files in SAM format as output. 
        --runThreads : Number Of Threads
        --runMode : genomeGenerates option which directs STAR to run genome indices generation job
        --genomeDir : specifies path to the directory
        --readFilesIn /path to read1 : indicates the path and the files containing the sequences to be mapped (i.e. RNA-seq FastQ files). 
        --outFileNamePrefix : by default, we have : ./ 
        The string indicates the output files name prefix (including pull or relative path).
        --outSAMtype : by default, we have : SAM ; the string indicates the type of SAM/BAM output
        BAM : 1st word : the output is a file in BAM format without sorting 
        SortedByCoordinate :  2nd word : this option allocates extra memory for sorting, which can be specified by -limitBAMtoRAM
        --outFilterMismatchNoverLmax : max number of mismatches relative to the read length
        if 75bp and we want 2 mismatch max : 0.027 
        if 50bp and we want 2 mismatch max : 0.04
        */

        """
        STAR --runThreadN ${cpus} \
        --runMode alignReads \
        --outFilterMismatchNoverLmax ${lmax} \
        --genomeDir ./index/ \
        --readFilesIn ${reads} \
        --outFileNamePrefix ./${file_id} \
        --outSAMtype BAM SortedByCoordinate 
        samtools index -@ ${cpus} -b ${file_id}*.out.bam ${file_id}_star.bam.bai 
        mv ${file_id}Log.final.out ${file_id}_star.statsLog.final.out
        """
    }
}
 
/* 
* Initialization of default empty values if the only the mapper STAR is used.
*/
if (params.onlySTAR || params.skipMapping){
  Channel
               .empty()
               .into { bowtie2_bam_files ; bowtie2_mapping_report}
}
/* 
* Alignment of reads using Bowtie2 mapper.
*/
else {
  process Bowtie2 {
      label "bowtie2"
          tag "$file_id"
          publishDir "${params.outdir}/mapping/BOWTIE2/$file_id", mode: 'copy'
            
          input:
          set file_id, file(reads) from fastq_trim_files2
          file index from index_files_BOWTIE.collect()

          output:
          set file_id, val(data_type), "*_sorted.bam" into bowtie2_bam_files
          file "*bowtie2.stats.log" into bowtie2_mapping_report

          script:
          data_type="_bowtie2_"
          index_id = index[0]
          for (index_file in index) {
            if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
              index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
            }
          }
          
          /*
           bowtie2 takes a Bowtie 2 index and a set of sequencing read files and outputs a set of alignments in SAM format.
           -P number of threads
           -x index file
           -U reads files
           -S output file in sam format
           -D number of seeding attempt that can fail before bowtie2 moves on another seed : 20 allow us to be very sensitive , but we can reduce execution time by lowering this number -> 15
           -R maximum number of times bowtie2 will try to "re-seed" a read with repetitive seeds ( for more information on reads with repetitive seeds , look for the manual of bowtie2). We use the default value here. 
           -N Number of missmatchs allowed in a seed alignement 
           -L length od the seeds substrings during the multiseed alignement : small value make alignement more sensitive but slower. default is 22 , we choose 20 to be more sensitive.
           -i function used to interval between seeds substrings during the multistring alignement. S,1,0.50 means : f(x)=1+0.5*sqrt(x) where x is the read size. this value allow for a lot of seeds substring to be used during multiseeds alignement. We can reduce the execution time (while bbeing less sensitive) bye choosing another function for exemple : S,1,1.15
          
          The grep command allows to select reads only aligned with 0, 1 ou 2 mismatches  

          Samtools allows to successively convert sam files into bam files, sort the bam files and index the sorted bam files
          */

          """
          bowtie2 -p ${cpus} \
          -D 20 -R 2 -N 1 -L 20 -i S,1,0.50 \
          -x ${index_id} \
          -U ${reads} \
          -S ${file_id}_bowtie2_tmp.sam 2> ${file_id}_bowtie2.stats.log
          
          grep -e "^@" -e 'XM:i:[012]' ${file_id}_bowtie2_tmp.sam > ${file_id}_bowtie2.sam

          samtools view -bS ${file_id}_bowtie2.sam -o ${file_id}_bowtie2.bam
          samtools sort -@ ${cpus} -O BAM -o ${file_id}_bowtie2_sorted.bam ${file_id}_bowtie2.bam
          samtools index -@ ${cpus} -b ${file_id}_bowtie2_sorted.bam ${file_id}_bowtie2_sorted.bam.bai
          """
  }
}

///////////////////////////////////////////////////////////////////////////////
/* --                         VARIANT CALLING                             -- */
///////////////////////////////////////////////////////////////////////////////

/* 
* Concat the sorted bam files from STAR and Bowtie2 alignment into a single 
* Channel is both mapper were used.
*/
if (!params.onlySTAR && !params.onlyBowtie2){
  bowtie2_bam_files.concat(star_bam_files).into { bam_files ; bam_files_to_bai }
}
/* 
* Set the sorted bam file to only the ones from the STAR alignment if this 
* option is chosen. 
*/
if (params.onlySTAR){
  star_bam_files.into { bam_files ; bam_files_to_bai }
}
/* 
* Set the sorted bam file to only the ones from the Bowtie2 alignment if this 
* option is chosen. 
*/
if (params.onlyBowtie2){
  bowtie2_bam_files.into { bam_files ; bam_files_to_bai }
}

/* 
* Variant Calling if mapping step is not skipped.
*/
if (!params.skipMapping){
  /*
  * Variant calling process that write the results into .cvs files
  */
  process Bcftools {
      label "bcftools"
          tag "$file_id"
          publishDir "${params.outdir}/bcftools/$file_id", mode: 'copy'
  
          input:
          set fasta_id, file(fasta) from fasta_2variantCalling.collect()
          set file_id, data_type, file(reads) from bam_files
          
          output:
          set file_id, data_type, "${file_id}${data_type}*" into variant_calling_file
          file "*.stats.txt" into bcftools_report 

          script:
          /*
          bcftools mpileup takes the aligments and the reference genome and generate a genotype likehoods at each genomic position
          -f is the reference genome in fasta format file 
          --no-BAQ; BAQ is the Phred-scaled probability of a read base being misaligned. 
          Applying this option greatly helps to reduce false SNPs caused by misalignments.

          bcftools call takes the genotype likehoods and generate the cause of the actual variance. 
          The output are goint to be the just the regions of the aligment which are different to our genome. 
          -m tells bcftool to use the default polling method. 
          -v output the variants sites only. -Ou output type u, which is the bcf format          
          
          bcftools stats Parse the bcf file and produce text file stats which is suitable for machine processing. 
          -F faidx indexed reference sequence file to determine INDEL context
          
          bcftools Query is a tools that allow you to view, subset and filter bcf file by quality. 
          -i include the variants, here with a threshole equal or more than 20. We only have variants that passed the filter. 
          It also re-format the vcf file into a text file for analysis, in this case a cvs file. 
          -H print header. 
          */
          """
          bcftools mpileup --no-BAQ -f ${fasta} ${reads} > ${file_id}${data_type}_pileup.vcf 
          bcftools call -mv -Ou -o ${file_id}${data_type}calls.bcf ${file_id}${data_type}_pileup.vcf 
          bcftools stats -F ${fasta} ${file_id}${data_type}calls.bcf > ${file_id}${data_type}calls.stats.txt 
          bcftools query -i 'QUAL>20 && DP>20' -f '%CHROM;%POS;%ID;%REF;%ALT;%QUAL;%FILTER\n' ${file_id}${data_type}calls.bcf -o ${file_id}${data_type}_calls.csv -H
          """
  }
}

///////////////////////////////////////////////////////////////////////////////
/* --                      MERGE ALL STEPS REPORTS                        -- */
///////////////////////////////////////////////////////////////////////////////

/*
* Regroups all reports from FastQC, Trimmomatic, STAR and Bowtie2 into a single
* report with MultiQC.
*/
process MultiQC {
    label "multiQC"
      publishDir "${params.outdir}/multiQC", mode: 'copy'

      input:
      file report_fastqc from fastqc_report.collect().ifEmpty([])
      file report_trim from trimming_report.collect().ifEmpty([])
      file report_star from star_mapping_report.collect().ifEmpty([])
      file report_mapping from bowtie2_mapping_report.collect().ifEmpty([])
      file report_bcftools from bcftools_report.collect().ifEmpty([])

      output:
      file "*multiqc_*" into multiqc_report

      when:
      !params.skipMultiqc

      script:
      """
      multiqc . 
      """
}
