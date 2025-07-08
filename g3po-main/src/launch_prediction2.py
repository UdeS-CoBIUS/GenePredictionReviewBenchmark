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

            ligne = ligne.strip().split(",")
            espece = ligne[0]
            print(ligne)
            if program == "augustus":
                code = ligne[1]

                if code not in dico_code.keys():
                    dico_code[code] = []
                    dico_code[code].append(espece)
                else:
                    dico_code[code].append(espece)
            elif program == "helixer":
                code = ligne[6] 
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
    # Configuration préalable
    # os.environ['AUGUSTUS_CONFIG_PATH'] ='/data/goglu/data/djoa2702/G3PO/g3po-main/src/augustus/config/'
    
    if add_nuc == 150:
        output = work_path + f"Predictions/augustus/150bp/"
        fasta_dir = work_path + f"References/Fasta_confirmed/150bp/"
    else:
        output = work_path + f"Predictions/augustus/{str(add_nuc)}Kb/"
        fasta_dir = work_path + f"References/Fasta_confirmed/{str(add_nuc)}Kb/"

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

        #break
    t1 = time.time()

    t_final = t1 - t0
    time_sortie.write(f"Augustus\t{add_nuc}\t{t_final}\n")
    print(f"Augustus took {t_final} seconds ")
    print(f"{add_nuc}Kb")

def launch_helixer(work_path, add_nuc, path_soft):
    """
    Launch Helixer for gene prediction.
    :param work_path: working folder path (root)
    :param add_nuc: Number of flanked nucleotides (not used directly in Helixer but kept for consistenc>
    :param path_soft: Directory where Helixer is installed
    :return:
    """
    # Define the output and fasta directories based on add_nuc parameter
    if add_nuc == 150:
        output = work_path + f"Predictions/helixer/150bp/"
        output_protein = work_path + f"Predictions/Proteines/helixer/150bp/"
        fasta_dir = work_path + f"References/Fasta_confirmed/150bp/"

    else:
        output = work_path + f"Predictions/helixer/{str(add_nuc)}Kb/"
        output_protein = work_path + f"Predictions/Proteines/helixer/{add_nuc}Kb/"
        fasta_dir = work_path + f"References/Fasta_confirmed/{str(add_nuc)}Kb/"
    
    # Ensure output directory exists
    os.makedirs(output, exist_ok=True)
    
    # Get the list of fasta files to process
    listf = os.listdir(fasta_dir)
    

    # Dictionary mapping species to their respective models
    dico_code = create_dico_model(work_path, "helixer")
    
    t0 = time.time()
    for elmt in tqdm(listf):
        nom = elmt.split("_")[1].replace(".fasta", "")
        for k, v in dico_code.items():
            if nom in v:
                species = k
                # Construct the command to run Helixer
                cmd = f"python {path_soft}/Helixer.py --lineage {species} --fasta-path {fasta_dir + elmt} --gff-output-path {output}helixer_{elmt}.gff3"
                os.system(cmd)
        #break
    
    t1 = time.time()
    time_sortie.write(f"Helixer\t{add_nuc}\t{t1 - t0}\n")
    print(f"Helixer took {t1 - t0} seconds for {add_nuc}Kb")

def launch_hmmgene(workpath, add_nuc, path_soft):
    work_path = os.getcwd() + "/"
    work_path = work_path.replace("src/", "")

    if add_nuc == 150:
        output = workpath + f"Predictions/hmmgene/150bp/"
        fasta_dir = workpath + f"References/Fasta_confirmed/150bp/"
    else:
        output = workpath + f"Predictions/hmmgene/{str(add_nuc)}Kb/"
        fasta_dir = workpath + f"References/Fasta_confirmed/{str(add_nuc)}Kb/"

    os.makedirs(output, exist_ok=True)
    listf = os.listdir(fasta_dir)
    #dico_code = create_dico_model(workpath, "geneid")
    t0 = time.time()

    for elmt in tqdm(listf):

        os.system(f"perl /usr/local/bin/hmmgene {fasta_dir}{elmt} > {output}hmmgene_{elmt}")
        #break
    os.chdir(work_path)

    t1 = time.time()

    t_final = t1 - t0
    time_sortie.write(f"HMMgene\t{add_nuc}\t{t_final}\n")
    print(f"HMMgene took {t_final} seconds ")
    print(f"{add_nuc}Kb")


if __name__ == '__main__':
    # workpath is the root directory where you download the git repository
    workpath =  os.getcwd() + "/../"

    # Please change the empty fields by the path of the different programs (eg: "/home/user/Genscanlinux/"")
    path_soft_augustus = "/usr/"
    path_soft_helixer = "/data/goglu/data/djoa2702/G3PO/g3po-main/src/helixer/Helixer/"
    path_soft_hmmgene = "/data/goglu/data/djoa2702/G3PO/g3po-main/src/hmmgene/"


    tt0 = time.perf_counter()
    _nuc = [150,2,4,6,8,10]
    time_sortie = open("./time_prediction.csv", "w")

    for add_nuc in _nuc:
        #launch_augustus(workpath, add_nuc, path_soft_augustus)
        #launch_helixer(workpath, add_nuc, path_soft_helixer) 
        launch_hmmgene(workpath, add_nuc, path_soft_hmmgene)

    tt1 = time.perf_counter()

    final_time = tt1-tt0
    
    time_sortie.close()
    print(f"time clock : {final_time}")
