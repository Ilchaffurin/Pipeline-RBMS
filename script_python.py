#! /usr/bin/python
#-*-coding: utf-8-*-


# commencer par taper la commande pour convertir le fichier au format BCF en un fichier au format VCF
# bcftools view -v output.bcf > output.vcf
 
#import sys

def conversion(f_entrance, f_exit_write):
    with open(f_entrance, "r") as file_vcf:
        dict_general = {}
        liste_generale = []
        for line in file_vcf:
            print(line)
            petit_dict_vcf = {}
            if line.startswith("##"):
                continue
            else:
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


    # Ecriture du fichier au format csv

    with open(f_exit_write, "w") as file_csv_write:
        for i in liste_generale:
            petit_dict_vcf = i
            file_csv_write.write("{0};{1};{2};{3};{4};{5};{6}\n".format(petit_dict_vcf["ID"], petit_dict_vcf["CHROM"], petit_dict_vcf["POS"], petit_dict_vcf["REF"], petit_dict_vcf["ALT"], petit_dict_vcf["QUAL"], petit_dict_vcf["FILTER"]))
        print(file_csv_write)
    print("Fin de l'écriture du fichier .csv")
print(conversion("calls.vcf", "calls.csv"))
        