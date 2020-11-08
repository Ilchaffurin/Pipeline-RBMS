#!/usr/bin/env nextflow


//////// OBJECTIF : Faire un pipeline pour identifier et analyser les variations génétiques des séquences nucléotiques des ARNr 

/////// Partie conversion du fichier VCF en fichier csv



///////// Commande : bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' calls.vcf -o file.csv -H
//////// ou : bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' calls.vcf -o file.csv
// Normalement, possibilité de parser le fichier de sortie (fichier au format csv) avec les headers, mais à voir... à tester lors du parsing du fichier



process Conversion {
    input:
    output:
    """

    """
}