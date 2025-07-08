#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created : 25/08/2019

@author: Scalzitti Nicolas 
"""

import os
from tqdm import tqdm

def extract_prot_seq_augustus(work_path, add_nuc):
    """
    From augustus prediction file, the function extract all predicted protein sequence and store them in output
    1 sequence = 1 file
    :param work_path: The root of the main project
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :return:
    """

    soft = "augustus"

    if add_nuc != 150:
        output = work_path + f"Predictions/Proteines/{soft}/{str(add_nuc)}Kb/"
        os.makedirs(output, exist_ok=True)
        path_pred = work_path + f"Predictions/{soft}/{str(add_nuc)}Kb/"

    else:
        add_nuc = int(add_nuc)
        output = work_path + f"Predictions/Proteines/{soft}/150bp/"
        os.makedirs(output, exist_ok=True)
        path_pred = work_path + f"Predictions/{soft}/150bp/"

    liste_protein_pred = os.listdir(path_pred)

    for protein in tqdm(liste_protein_pred):

        borne = []
        liste_borne = []

        with open(path_pred + protein, "r") as file_R1:
            for i, ligne in enumerate(file_R1):
                if "[" in ligne:
                    a = i
                    borne.append(a)
                if "]" in ligne:
                    b = i
                    borne.append(b)

                if len(borne) == 2:
                    liste_borne.append(borne)
                    borne = []

        proteins = []
        with open(path_pred + protein, "r") as file_R2:
            for i, ligne in enumerate(file_R2):
                for bornes in liste_borne:
                    for j in range(bornes[0], bornes[1] + 1):
                        if i == j:
                            proteins.append(ligne)

        num = 1

        for p in proteins:
            # 1st line
            if "[" in p:
                pro = []
                pro.append(p)
            # last line
            elif "]" in p:
                pro.append(p)
                x = ",".join(pro)

                new_prot = str(x).replace("\n", "").replace("# protein sequence = [", "").replace(",# ", "").replace("]", "")
                with open(output + f"{protein.replace('.fasta', '')}_{str(num)}.fasta", "w") as file_W1:
                    file_W1.write(new_prot)

                pro = []
                # Increment the number of the predicted protein
                num += 1

            # Intermediate line
            else:
                pro.append(p)

        for p in proteins:
            if "[" in p and "]" in p:
                new_prot = str(p).replace("\n", "").replace("# protein sequence = [", "").replace(",# ", "").replace("]", "")
                with open(output + f"{protein.replace('.fasta', '')}_{str(num)}.fasta", "w") as file_W1:
                    file_W1.write(new_prot)

def extract_prot_seq_genscan(work_path, add_nuc):
    """
    From genscan prediction file, the function extract all predicted protein sequence and store them in output
    1 sequence = 1 file
    :param work_path: The root of the main project
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :return:
    """
    soft = "genscan"

    if add_nuc != 150:
        output = work_path + f"Predictions/Proteines/{soft}/{str(add_nuc)}Kb/"
        os.makedirs(output, exist_ok=True)
        path_pred = work_path + f"Predictions/{soft}/{str(add_nuc)}Kb/"

    else:
        add_nuc = int(add_nuc)
        output = work_path + f"Predictions/Proteines/{soft}/150bp/"
        os.makedirs(output, exist_ok=True)
        path_pred = work_path + f"Predictions/{soft}/150bp/"

    liste_protein_pred = os.listdir(path_pred)

    for protein in tqdm(liste_protein_pred):
        with open(path_pred + protein, "r") as file_R1:
            for i, ligne in enumerate(file_R1):

                if ">" in ligne:
                    a = file_R1.readlines()

                if "NO EXONS/GENES PREDICTED IN SEQUENCE" in ligne:
                    a = []

        a.append("\n")
        num = 1
        pro = []

        for p in a:
            if "\n" == p:
                x = ",".join(pro)
                with open(output + protein.replace(".fasta", "") + "_" + str(num) + ".fasta", "w") as file_W1:
                    file_W1.write(str(x.replace("\n", "").replace(",", "")))
                pro = []
                num += 1

            elif ">" in p:
                pass
            else:
                pro.append(p)

def extract_prot_seq_geneid(work_path, add_nuc):
    """
    From geneid prediction file, the function extract all predicted protein sequence and store them in output
    1 sequence = 1 file
    :param work_path: The root of the main project
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :return:
    """
    soft = "geneid"

    if add_nuc != 150:
        output = work_path + f"Predictions/Proteines/{soft}/{str(add_nuc)}Kb/"
        os.makedirs(output, exist_ok=True)
        path_pred = work_path + f"Predictions/{soft}/{str(add_nuc)}Kb/"

    else:
        add_nuc = int(add_nuc)
        output = work_path + f"Predictions/Proteines/{soft}/150bp/"
        os.makedirs(output, exist_ok=True)
        path_pred = work_path + f"Predictions/{soft}/150bp/"

    liste_protein_pred = os.listdir(path_pred)

    for protein in tqdm(liste_protein_pred):
        proteine = []

        with open(path_pred + protein, "r") as file_R1:
            for i, ligne in enumerate(file_R1):
                if "#" in ligne or "\t" in ligne:
                    pass
                else:
                    proteine.append(ligne)
        num = 1
        pro = []
        for p in proteine:

            if "\n" == p:
                x = ",".join(pro)
                if len(x) > 0:
                    with open(output + protein.replace(".fasta", "") + "_" + str(num) + ".fasta", "w") as file_W1:
                        file_W1.write(str(x.replace("\n", "").replace(",", "")))
                    pro = []
                    num += 1
            elif ">" in p:
                pass
            else:
                pro.append(p)

def extract_prot_seq_glimmer(work_path, add_nuc):
    """
    From glimmer prediction file, the function extract all predicted protein sequence and store them in output
    all sequences = 1 file. The glimmer_2_aa.py script must be placed in the same directory.
    :param work_path: The root of the main project
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :return:
    """
    soft  = "glimmer"
    path_script = os.getcwd()

    if add_nuc != 150:
        liste_proteine = os.listdir(work_path + f"References/Fasta/{str(add_nuc)}Kb/")
        output = work_path + f"Predictions/Proteines/{soft}/{str(add_nuc)}Kb/"
        os.makedirs(output, exist_ok=True)
    else:
        liste_proteine = os.listdir(work_path + f"References/Fasta/150bp/")
        add_nuc = int(add_nuc)
        output = work_path + f"Predictions/Proteines/{soft}/150bp/"
        os.makedirs(output, exist_ok=True)

    for proteine in tqdm(liste_proteine):
        proteine, ext = os.path.splitext(proteine)

        if add_nuc != 150:
            predictions = work_path + f"Predictions/glimmer/{str(add_nuc)}Kb/glimmer_{proteine}.fasta"
            fasta = work_path + f"References/Fasta/{str(add_nuc)}Kb/{proteine}"
        else:
            predictions = work_path + f"Predictions/glimmer/150bp/glimmer_{proteine}.fasta"
            fasta = work_path + f"References/Fasta/150bp/{proteine}"

        os.system(f"python2 {path_script}/glimmer_2_aa.py {work_path} {str(add_nuc)} {predictions} {fasta}.fasta")



if __name__ == '__main__':
    # workpath is the root directory where you download the git repository
    workpath = "my_path/Benchmark_study/"
    _nuc = [150,2,4,6,8,10]

    for add_nuc in _nuc:
        extract_prot_seq_augustus(workpath, add_nuc)
        extract_prot_seq_genscan(workpath, add_nuc)
        extract_prot_seq_geneid(workpath, add_nuc)
        extract_prot_seq_glimmer(workpath, add_nuc)




