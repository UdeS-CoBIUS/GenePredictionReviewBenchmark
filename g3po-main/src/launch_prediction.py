#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created : 20/08/2019

@author: Scalzitti Nicolas 
"""

import os
from tqdm import tqdm
import time

def create_dico_model(work_path, program):
    """
    Generate a dictionnary with the ID_specie and the model associate, Key = code, value = [spe1, spe2...]
    :param work_path: working folder path (root)
    :param program: Name of the program
    :return:
    """
    dico_code = {}

    with open(work_path + f"Models/models.csv", "r") as file_R1:
        for i, ligne in enumerate(file_R1):

            ligne = ligne.strip().split(";")
            espece = ligne[0]

            if program == "augustus":
                code = ligne[1]

                if code not in dico_code.keys():
                    dico_code[code] = []
                    dico_code[code].append(espece)
                else:
                    dico_code[code].append(espece)
            elif program == "genscan":
                code = ligne[2]

                if code not in dico_code.keys():
                    dico_code[code] = []
                    dico_code[code].append(espece)
                else:
                    dico_code[code].append(espece)
            elif program == "geneid":
                code = ligne[3]

                if code not in dico_code.keys():
                    dico_code[code] = []
                    dico_code[code].append(espece)
                else:
                    dico_code[code].append(espece)
            elif program == "glimmer":
                code = ligne[4]

                if code not in dico_code.keys():
                    dico_code[code] = []
                    dico_code[code].append(espece)
                else:
                    dico_code[code].append(espece)
            elif program == "snap":
                code = ligne[5]

                if code not in dico_code.keys():
                    dico_code[code] = []
                    dico_code[code].append(espece)
                else:
                    dico_code[code].append(espece)

    return dico_code


def launch_augustus(work_path, add_nuc, path_soft):
    """
    Launch the soft Augustus
    :param work_path: working folder path (root)
    :param add_nuc: Number of flanked nucleotids
    :param path_soft : Directory where are installed Augustus
    :return:
    """

    if add_nuc == 150:
        output = work_path + f"Predictions/augustus/150bp/"
        fasta_dir = work_path + f"References/Fasta/150bp/"
    else:
        output = work_path + f"Predictions/augustus/{str(add_nuc)}Kb/"
        fasta_dir = work_path + f"References/Fasta/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)

    # Directory of the .fasta files
    listf = os.listdir(fasta_dir)
    dico_code = create_dico_model(work_path, "augustus")

    t0 = time.time()
    for elmt in tqdm(listf):

        nom = elmt.split("_")[1].replace(".fasta", "")
        for k, v in dico_code.items():
            if nom in v:
                species = k
                os.system(f"{path_soft}bin/augustus --species={species} --softmasking=1 --gff3=off {fasta_dir + str(elmt)} > {output}augustus_{elmt}")

    t1 = time.time()

    t_final = t1 - t0
    print(f"Augustus took {t_final} seconds ")
    print(f"{add_nuc}Kb")


def launch_genscan(work_path, add_nuc, path_soft):
    """
    Launch the soft Genscan
    :param work_path: working folder path (root)
    :param add_nuc: Number of flanked nucleotids
    :param path_soft : Directory where are installed Genscan
    :return:
    """
    if add_nuc == 150:
        output = work_path + f"Predictions/genscan/150bp/"
        fasta_dir = work_path + f"References/Fasta/150bp/"
    else:
        output = work_path + f"Predictions/genscan/{str(add_nuc)}Kb/"
        fasta_dir = work_path + f"References/Fasta/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)

    listf = os.listdir(fasta_dir)

    dico_code = create_dico_model(work_path, "genscan")

    t0 = time.time()
    for elmt in tqdm(listf):
        nom = elmt.split("_")[1].replace(".fasta", "")
        for k, v in dico_code.items():
            if nom in v:
                species = k
                os.system(f" {path_soft}genscan {path_soft}{species} {fasta_dir}{elmt} > {output}genscan_{elmt}")

    t1 = time.time()

    t_final = t1 - t0
    print(f"Genscan took {t_final} seconds ")
    print(f"{add_nuc}Kb")


def launch_geneid(workpath, add_nuc, path_soft):
    """
    Launch the soft GeneID
    :param workpath: working folder path (root)
    :param add_nuc: Number of flanked nucleotids
    :param path_soft : Directory where are installed GeneID
    :return:
    """

    work_path = os.getcwd() + "/"
    work_path = work_path.replace("src/", "")

    if add_nuc == 150:
        output = workpath + f"Predictions/geneid/150bp/"
        fasta_dir = workpath + f"References/Fasta/150bp/"
    else:
        output = workpath + f"Predictions/geneid/{str(add_nuc)}Kb/"
        fasta_dir = workpath + f"References/Fasta/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)
    listf = os.listdir(fasta_dir)
    dico_code = create_dico_model(workpath, "geneid")
    t0 = time.time()

    for elmt in tqdm(listf):
        nom = elmt.split("_")[1].replace(".fasta", "")

        for k, v in dico_code.items():
            if nom in v:
                species = k
                os.chdir(path_soft)
                # New_param is a directory with new models for geneid 1.4 - Downloaded from their website and adding in
                # a new directory called 'new_param'
                os.system(f"bin/geneid -A -P param/new_param/{species} {fasta_dir}{elmt} > {output}geneid_{elmt}")
    os.chdir(work_path)

    t1 = time.time()

    t_final = t1 - t0
    print(f"Geneid took {t_final} seconds ")
    print(f"{add_nuc}Kb")


def launch_glimmer(work_path, add_nuc, path_soft):
    """
    Launch the soft Glimmer
    :param work_path: working folder path (root)
    :param add_nuc: Number of flanked nucleotids
    :param path_soft : Directory where are installed Glimmer
    :return:
    """
    if add_nuc == 150:
        output = work_path + f"Predictions/glimmer/150bp/"
        fasta_dir = work_path + f"References/Fasta/150bp/"

    else:
        output = work_path + f"Predictions/glimmer/{str(add_nuc)}Kb/"
        fasta_dir = work_path + f"References/Fasta/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)
    listf = os.listdir(fasta_dir)
    dico_code = create_dico_model(work_path, "glimmer")

    t0 = time.time()
    for elmt in tqdm(listf):
        nom = elmt.split("_")[1].replace(".fasta", "")
        for k, v in dico_code.items():
            if nom in v:
                species = k
                os.system(
                    f"{path_soft}bin/glimmerhmm_linux_x86_64 {fasta_dir}{elmt} -d {path_soft}trained_dir/{species} -g > {output}glimmer_{elmt}")
    t1 = time.time()

    t_final = t1 - t0
    print(f"Glimmer took {t_final} seconds ")
    print(f"{add_nuc}Kb")


def launch_snap(work_path, add_nuc, path_soft):
    """
    Launch the soft Snap, predicted proteins are stored in /Predictions/Proteines/snap/{add_nuc}Kb/
    :param work_path: working folder path (root)
    :param add_nuc: Number of flanked nucleotids
    :param path_soft : Directory where are installed Snap
    :return:
    """
    if add_nuc == 150:
        output = work_path + f"Predictions/snap/150bp/"
        output_protein = work_path + f"Predictions/Proteines/snap/150bp/"
        fasta_dir = work_path + f"References/Fasta/150bp/"

    else:
        output = work_path + f"Predictions/snap/{str(add_nuc)}Kb/"
        output_protein = work_path + f"Predictions/Proteines/snap/{add_nuc}Kb/"
        fasta_dir = work_path + f"References/Fasta/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)
    os.makedirs(output_protein, exist_ok=True)
    listf = os.listdir(fasta_dir)
    dico_code = create_dico_model(work_path, "snap")

    t0 = time.time()
    for elmt in tqdm(listf):
        nom = elmt.split("_")[1].replace(".fasta", "")
        for k, v in dico_code.items():
            if nom in v:
                species = k
                os.system(
                    f"{path_soft}snap -gff -quiet -lcmask {path_soft}HMM/{species} {fasta_dir}{elmt} -aa {output_protein}{elmt} > {output}snap_{elmt}")
    t1 = time.time()

    t_final = t1 - t0
    print(f"Snap took {t_final} seconds ")
    print(f"{add_nuc}Kb")


if __name__ == '__main__':
    # workpath is the root directory where you download the git repository
    workpath =  os.getcwd() + "/../"

    # Please change the empty fields by the path of the different programs (eg: "/home/user/Genscanlinux/"")
    path_soft_augustus = ""
    path_soft_geneid = ""
    path_soft_genscan = ""
    path_soft_glimmer = ""
    path_soft_snap = ""

    tt0 = time.clock()
    _nuc = [150,2,4,6,8,10]

    for add_nuc in _nuc:
        launch_augustus(workpath, add_nuc, path_soft_augustus)
        launch_geneid(workpath, add_nuc, path_soft_geneid)
        launch_genscan(workpath, add_nuc, path_soft_genscan)
        launch_glimmer(workpath, add_nuc, path_soft_glimmer)
        launch_snap(workpath, add_nuc, path_soft_snap)

    tt1 = time.clock()

    final_time = tt1-tt0
    print(f"time clock : {final_time}")