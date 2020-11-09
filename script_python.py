#!/usr/bin/python
#-*-coding: utf-8-*-

# Conversion file in VFC format into file in csv format
# VCF format looks like BED one

import numpy as np
import vcf
import csv

def conversion(f_entrance, f_exit):
    with open(f_entrance, "r") as file_vcf:
        file_vcf = f_entrance.read()
        #print(file_vcf)
        liste_generale = []
        for line in file_vcf:
            petit_dict_vcf = {}
            colonne = line.split("\t")
            petit_dict_vcf["CHROM"] = colonne[0]
            petit_dict_vcf["POS"] = colonne[1]
            petit_dict_vcf["ID"] = colonne[2]
            petit_dict_vcf["REF"] = colonne[3]
            petit_dict_vcf["ALT"] = colonne[4]
            petit_dict_vcf["QUAL"] = colonne[5]
            petit_dict_vcf["FILTER"] = colonne[6]
            # on ajoute le petit dictionnaire dans la liste générale
            liste_generale.append(petit_dict_vcf)
        print(liste_generale)
    print("Fin de la lecture du fichier .vcf")





    # Writing to file in csv format

    with open(f_exit, "w", nexline = '') as file_csv:
        writerCSV = csv.writer(file_csv, delimiter = "\t")
        writerCSV = writerow(["[1]CHROM", "[2]POS", "[3]ID", "[4]REF", "[5]ALT", "[6]QUAL", "[7]FILTER"])