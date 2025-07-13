#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created : 27/08/2019

@author: Scalzitti Nicolas 
"""

import os
from tqdm import tqdm

def moyenne(tableau):
    return sum(tableau) / len(tableau)

def liste_gene():
    liste_validated = []

    with open("/home/scalzitti/These/G3PO_Benchmark/Benchmark_study/Sequence_types/All.seq", "r") as file_R0:
        for ligne in file_R0:
            ligne = ligne.strip()
            liste_validated.append(ligne)
    return liste_validated


def join_augustus_genscan_geneid(work_path, soft, add_nuc):
    """
    Join reference and predicted protein sequence for augustus, genscan and geneid program
    :param work_path: The root of the main project
    :param soft: Name of the predictor
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :return:
    """

    if add_nuc == 150:
        path_pred = work_path + f"Predictions/Proteines/{soft}/150bp/"
        path_ref = work_path + "References/Proteins/"
        output = work_path + f"Predictions/Proteines/Joined/{soft}/150bp/"
    else:
        path_pred = work_path + f"Predictions/Proteines/{soft}/{str(add_nuc)}Kb/"
        path_ref = work_path + "References/Proteins/"
        output = work_path + f"Predictions/Proteines/Joined/{soft}/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)

    liste_ref = os.listdir(path_ref)
    liste_pred = os.listdir(path_pred)

    for file_ref in tqdm(liste_ref):
        for file_pred in liste_pred:

            filename = file_ref.split(".")[0]

            if filename in file_pred:
                with open(output + file_pred, "a") as file_W1:
                    file_W1.write("> REFERENCE " + file_ref + "\n")
                    with open(path_ref + file_ref, 'r') as file_R1:
                        for ligne in file_R1:
                            if ">" in ligne:
                                pass
                            else:
                                l1 = (ligne.replace("\n", ""))
                                file_W1.write(l1)

                    file_W1.write("\n> PREDICTION " + file_pred + "\n")
                    with open(path_pred + file_pred, 'r') as file_R2:
                        for ligne in file_R2:
                            l2 = ligne
                            file_W1.write(l2)


def join_glimmer_snap(work_path, soft, add_nuc):
    """
    Join reference and predicted protein sequence for glimmer and snap program
    :param work_path: The root of the main project
    :param soft: Name of the predictor
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :return:
    """

    if add_nuc == 150:
        path_pred = work_path + f"Predictions/Proteines/{soft}/150bp/"
        path_ref = work_path + "References/Proteins/"
        output = work_path + f"Predictions/Proteines/Joined/{soft}/150bp/"
    else:
        path_pred = work_path + f"Predictions/Proteines/{soft}/{str(add_nuc)}Kb/"
        path_ref = work_path + "References/Proteins/"
        output = work_path + f"Predictions/Proteines/Joined/{soft}/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)

    liste_proteine_ref = os.listdir(path_ref)
    liste_protein_pred = os.listdir(path_pred)

    # Write prediction
    for protein in tqdm(liste_protein_pred):

        seq_proteine = []
        prot = protein.replace(f"{soft}_", "").replace(".fasta", "")

        for file_ref in liste_proteine_ref:
            num = 0

            # Search the same ID
            if prot in file_ref:

                reference = []
                with open(path_ref + file_ref, "r") as file_r:
                    for ligne_r in file_r:
                        if ">" in ligne_r:
                            head = ligne_r
                            pass
                        else:
                            reference.append(ligne_r)

                # Protein sequence of the reference
                reference = ','.join(reference).replace(",", "").replace("\n", "")

                # last line of the prediction file
                with open(path_pred + protein, "r") as file_R1:
                    compteur = 0
                    for i in file_R1:
                        compteur += 1

                with open(path_pred + protein, "r") as file_R1:

                    for i, ligne in enumerate(file_R1):

                        if ">" in ligne:
                            head_pred = ligne
                            seq_proteine = ','.join(seq_proteine).replace(",", "").replace("\n", "")
                            if len(seq_proteine) == 0:
                                pass
                            else:
                                num += 1
                                with open(output + f"{prot}_{str(num)}.fasta", "w") as file_W1:
                                    file_W1.write(">REF_" + head.replace(">",""))
                                    file_W1.write(reference + "\n")
                                    file_W1.write(">PRED_" + head_pred.replace(">",""))
                                    file_W1.write(seq_proteine)
                            seq_proteine = []
                        else:
                            seq_proteine.append(ligne)
                            if i == compteur - 1:
                                num += 1
                                seq_proteine = ','.join(seq_proteine).replace(",", "").replace("\n", "")
                                with open(output + prot + "_" + str(num) + ".fasta", "w") as file_W1:
                                    file_W1.write(">REF_" + head.replace(">",""))
                                    file_W1.write(reference + "\n")
                                    file_W1.write(">PRED_" + head_pred.replace(">",""))
                                    file_W1.write(seq_proteine)


def launch_mafft(work_path, soft, add_nuc, mafft_path):
    """
    Launch the alignement with Mafft
    :param work_path: The root of the main project
    :param soft: Name of the predictor
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :param mafft_path: Mafft program path
    :return:
    """
    if add_nuc == 150:
        output = work_path + f"Predictions/Proteines/Alignements/{soft}/150bp/"
        path_couple = work_path + f"Predictions/Proteines/Joined/{soft}/150bp/"
    else:
        output = work_path + f"Predictions/Proteines/Alignements/{soft}/{str(add_nuc)}Kb/"
        path_couple = work_path + f"Predictions/Proteines/Joined/{soft}/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)
    liste_file = os.listdir(path_couple)

    for filename in tqdm(liste_file):
        print(">>> " + filename)
        os.system(f"{mafft_path} --auto --clustalout --reorder {path_couple}{filename} > {output}{filename}")


def similarity_augustus_genscan_geneid(work_path, soft, add_nuc):
    """
    Calculate the % similarity between reference and predicted protein sequence for the augustus/genscan/geneid
    alignement files
    :param work_path: the root of the main project
    :param soft: Name of the predictor
    :param add_nuc: number of nucleotides flanked upstream and downstream of gene
    :return:
    """

    if add_nuc == 150:
        path_res = work_path + f"Predictions/Proteines/Alignements/{soft}/150bp/"
        output = work_path + f"Resultats/Proteins/150bp/"
    else:
        path_res = work_path + f"Predictions/Proteines/Alignements/{soft}/{str(add_nuc)}Kb/"
        output = work_path + f"Resultats/Proteins/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)
    valeur_simil = []

    # File list with all gene name
    liste_file = os.listdir(path_res)

    with open(f"{output}{soft}_%similarity.csv", "w") as file_W1:
        # writing header
        file_W1.write("File;%_similarity;Len_prot\n")

        for filename in tqdm(liste_file):

            longueur_seq = 0
            long_seq_pred = 0
            nombre_star = 0

            with open(path_res + filename, "r") as file_R1:
                for ligne in file_R1:
                    # add sequence length
                    if "REFERENCE" in ligne:
                        ligne = ligne.strip().split(" ")
                        ligne = list(filter(None, ligne))
                        longueur_seq += len(ligne[1])

                    if "PREDICTION" in ligne:
                        ligne = ligne.strip().split(" ")
                        ligne = list(filter(None, ligne))
                        for lettre in ligne[1]:
                            if lettre == "-":
                                pass
                            else:
                                long_seq_pred += 1

                    if "REFERENCE" in ligne or "PREDICTION" in ligne:
                        pass
                    else:
                        for lettre in ligne:
                            if lettre == "*":
                                nombre_star += 1

            similarite = (nombre_star * 100) / longueur_seq

            file_W1.write(f"{filename};{str(round(similarite, 2))};{str(longueur_seq)}\n")
            valeur_simil.append(round(similarite, 2))

        m2 = moyenne(valeur_simil)
        print(f"Mean for all: {m2}")

    l_f = liste_gene()
    l_gene_avec_simil = set()

    with open(f"{output}{soft}_%similarity.csv", "r") as file_R3:
        for i, ligne in enumerate(file_R3):
            if i != 0:
                id_prot = ligne.strip().split(";")
                id_prot = id_prot[0].split("_")
                id_prot = id_prot[1] + "_" + id_prot[2]
                l_gene_avec_simil.add(id_prot)

    for gene in l_f:
        if gene not in l_gene_avec_simil:

            with open(f"{output}{soft}_%similarity.csv", "a") as file_W1:
                file_W1.write(f"{soft}_{gene}_0.fasta;0.0;0\n")

def similarity_glimmer_snap(work_path, soft, add_nuc):
    """
    Calculate the % similarity between reference and predicted protein sequence for the glimmer/snap alignement files
    :param work_path: the root of the main project
    :param soft: Name of the predictor
    :param add_nuc: number of nucleotides flanked upstream and downstream of gene
    :return:
    """

    if add_nuc == 150:
        path_res = work_path + f"Predictions/Proteines/Alignements/{soft}/150bp/"
        output = work_path + f"Resultats/Proteins/150bp/"
    else:
        path_res = work_path + f"Predictions/Proteines/Alignements/{soft}/{str(add_nuc)}Kb/"
        output = work_path + f"Resultats/Proteins/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)
    valeur_simil = []

    liste_file = os.listdir(path_res)

    with open(f"{output}{soft}_%similarity.csv", "w") as file_W1:
        file_W1.write("File;%_similarity;Len_prot\n")
        for filename in tqdm(liste_file):

            longueur_seq = 0
            long_seq_pred = 0
            nombre_star = 0
            with open(path_res + filename, "r") as file_R1:
                for i, ligne in enumerate(file_R1):
                    if i == 0 or ligne == "\n":
                        pass
                    else:
                        if "REF" in ligne or "TrEMBL|" in ligne:

                            line = ligne.strip().split(" ")
                            line = list(filter(None, line))
                            longueur_seq += len(line[1])

                        name = filename.split("_")
                        name = name[1] + "_" + name[2]
                        if "PRED" in ligne:

                            ligne = ligne.strip().split(" ")
                            ligne = list(filter(None, ligne))
                            for lettre in ligne[1]:
                                if lettre == "-":
                                    pass
                                else:
                                    long_seq_pred += 1
                        else:
                            pass

                        for lettre in ligne:
                            if lettre == "*":
                                nombre_star += 1
            try:
                similarite = (nombre_star * 100) / longueur_seq
            except ZeroDivisionError:
                pass

            file_W1.write(f"{filename};{str(round(similarite, 2))};{str(longueur_seq)}\n")
            valeur_simil.append(round(similarite, 2))

        m2 = moyenne(valeur_simil)
        print(f"Mean for all: {m2}")
    l_f = liste_gene()

    l_gene_avec_simil = set()

    with open(f"{output}{soft}_%similarity.csv", "r") as file_R3:
        for i, ligne in enumerate(file_R3):
            if i != 0:
                id_prot = ligne.strip().split(";")
                id_prot = id_prot[0].split("_")
                id_prot = id_prot[0] + "_" + id_prot[1]
                l_gene_avec_simil.add(id_prot)

    for gene in l_f:
        if gene + ".fasta" not in l_gene_avec_simil:
            with open(f"{output}{soft}_%similarity.csv", "a") as file_W1:
                file_W1.write(f"{gene}_0.fasta;0.0;0\n")


def best_simmilarity(work_path, soft, add_nuc):
    """
    Calculate the best % similarity between reference and predicted protein sequence.
    :param work_path: the root of the main project
    :param soft: Name of the predictor
    :param add_nuc: number of nucleotides flanked upstream and downstream of gene
    :return:
    """

    if add_nuc == 150:
        path = work_path + f"Resultats/Proteins/150bp/"
        path_res = f"{path}{soft}_%similarity.csv"
        output = work_path + f"Resultats/Proteins/150bp/"
    else:
        path = work_path + f"Resultats/Proteins/{str(add_nuc)}Kb/"
        path_res = f"{path}{soft}_%similarity.csv"
        output = work_path + f"Resultats/Proteins/{str(add_nuc)}Kb/"

    with open(path_res, "r") as file_R1:
        dico_filename = {}
        for i, ligne in enumerate(file_R1):
            if i != 0:
                try:
                    filename = ligne.strip().split(";")[0]
                    filename = filename.split("_")
                    len_prot = ligne.strip().split(";")[2]

                    if soft == 'glimmer' or soft == 'snap':
                        num = filename[2].split(".")
                    else:
                        num = filename[3].split(".")

                    num = int(num[0])

                    if soft == 'glimmer' or soft == 'snap':
                        filename = filename[0] + "_" + filename[1]
                    else:
                        filename = filename[1] + "_" + filename[2]

                    similarity = float(ligne.strip().split(";")[1])

                    if filename not in dico_filename.keys():
                        dico_filename[filename] = []
                        dico_filename[filename].append([similarity, num, len_prot])

                    else:
                        dico_filename[filename].append([similarity, num, len_prot])

                except IndexError:
                    pass

        mean_all = []
        nbr_100 = 0

        with open(f"{output}{soft}_best_%similarity.csv", "w") as file_W1:
            for k, v in dico_filename.items():
                num = max(v)[1]
                similarity = max(v)[0]
                len_prot = max(v)[2]

                if similarity == 100.0:
                    nbr_100 += 1

                mean_all.append(similarity)
                file_W1.write(f"{k + '_' + str(num)};{v};{similarity};{int(len_prot)}\n")
            file_W1.write(f"{round(moyenne(mean_all), 3)};{nbr_100}")


def main(workpath, mafft_path):
    l_soft = ["augustus", "genscan", "geneid", "glimmer", "snap"]
    _nuc = [150, 2, 4, 6, 8, 10]
    l_soft_agg = ["augustus", "genscan", "geneid"]
    l_soft_gs = ["glimmer", "snap"]



    print(">>> protein sequences joined")
    for soft in l_soft_agg:
        for add_nuc in _nuc:
            join_augustus_genscan_geneid(workpath, soft, add_nuc)

    for soft in l_soft_gs:
        for add_nuc in _nuc:
            join_glimmer_snap(workpath, soft, add_nuc)

    print(">>> Mafft Launch")
    for soft in l_soft:
        for add_nuc in _nuc:
            launch_mafft(workpath, soft, add_nuc, mafft_path)

    print(">>> Calculation of identity percentage")
    for soft in l_soft_agg:
        for add_nuc in _nuc:
            similarity_augustus_genscan_geneid(workpath,soft, add_nuc)

    for soft in l_soft_gs:
        for add_nuc in _nuc:
            similarity_glimmer_snap(workpath,soft, add_nuc)

    for soft in l_soft:
        for add_nuc in _nuc:
            best_simmilarity(workpath, soft,add_nuc)

if __name__ == '__main__':
    # workpath is the root directory where you download the git repository
    workpath = "my_path/Benchmark_study/"
    mafft_path = "mafft_path"

    main(workpath, mafft_path)




