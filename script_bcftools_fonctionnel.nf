#!/usr/bin/env nextflow


process Conversion {
    label "conversion"
        tag "$file_id"
        publishDir ${calls.bcf}

        input:
        set file_id from variant_calling_file

        output:
        set file_id into csv_report

        script:
        """
        bcftools query -f '%CHROM;%POS;%ID;%REF;%ALT;%QUAL;%FILTER\n' calls.bcf -o file.csv -H
        """
}