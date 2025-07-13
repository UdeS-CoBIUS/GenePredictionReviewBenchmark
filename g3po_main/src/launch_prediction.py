#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created : 20/08/2019

@author: Scalzitti Nicolas 
@modified by Davy Ouedraogo on 09/07/2025
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
    path_to_models = os.path.join(work_path, "Models", "models.csv")

    with open(path_to_models, "r") as file_R1:
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

def launch_augustus(work_path, add_nuc, path_soft, test):
    """
    Launch the soft Augustus
    :param work_path: working folder path (root)
    :param add_nuc: Number of flanked nucleotids
    :param path_soft : Directory where are installed Augustus
    :return:
    """
    # Configuration préalable
    
    if test:
        try:
            output = os.path.join(os.getcwd(), "Test", "Predictions", "augustus", "150bp")
            os.makedirs(output, exist_ok=True)
            #os.system("/usr/bin/augustus --species=leishmania_tarentolae --softmasking=1 --gff3=off /home/local/USHERBROOKE/ouew2201/Documents/GenePredictionReviewBenchmark/g3po_main/References/Fasta_confirmed/150bp/A4I4C4_LEIIN.fasta > /home/local/USHERBROOKE/ouew2201/Documents/GenePredictionReviewBenchmark/Test/Predictions/augustus/150bp/augustus_A4I4C4_LEIIN.fasta")
            os.system("/usr/bin/augustus --species=fly --softmasking=1 --gff3=off /home/local/USHERBROOKE/ouew2201/Documents/GenePredictionReviewBenchmark/g3po_main/References/Fasta_confirmed/150bp/B3M138_DROAN.fasta > /home/local/USHERBROOKE/ouew2201/Documents/GenePredictionReviewBenchmark/Test/Predictions/augustus/150bp/augustus_B3M138_DROAN.fasta")

                
            return True
        except:
            return False
    else:
        try:
            if add_nuc == 150:
                output = os.path.join(work_path, "Predictions", "augustus", "150bp")
                fasta_dir = os.path.join(work_path, "References", "Fasta_confirmed", "150bp")
            else:
                output = os.path.join(work_path, "Predictions", "augustus", f"{str(add_nuc)}Kb")
                fasta_dir = os.path.join(work_path, "References", "Fasta_confirmed", f"{str(add_nuc)}Kb")

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
                        os.system(f"{os.path.join(path_soft, 'augustus')} --species={species} --softmasking=1 --gff3=off {os.path.join(fasta_dir, str(elmt))} > {os.path.join(output, f'augustus_{elmt}')}")
                        #print(f"{os.path.join(path_soft, 'augustus')} --species={species} --softmasking=1 --gff3=off {os.path.join(fasta_dir, str(elmt))} > {os.path.join(output, f'augustus_{elmt}')}")
                #break
            t1 = time.time()

            t_final = t1 - t0
            #time_sortie.write(f"Augustus\t{add_nuc}\t{t_final}\n")
            print(f"Augustus took {t_final} seconds ")
            print(f"{add_nuc}Kb")
            return True
        except:
            return False

def launch_helixer(work_path, add_nuc, path_soft, test):
    """
    Launch Helixer for gene prediction.
    :param work_path: working folder path (root)
    :param add_nuc: Number of flanked nucleotides (not used directly in Helixer but kept for consistenc>
    :param path_soft: Directory where Helixer is installed
    :return:
    """
    # Define the output and fasta directories based on add_nuc parameter
    if test:
        try:
            output = os.path.join(os.getcwd(), "Test", "Predictions", "helixer", "150bp")
            os.makedirs(output, exist_ok=True)
            #os.system("python3.10 Helixer/Helixer.py --lineage invertebrate --fasta-path /home/local/USHERBROOKE/ouew2201/Documents/GenePredictionReviewBenchmark/g3po_main/References/Fasta_confirmed/150bp/A4I4C4_LEIIN.fasta --gff-output-path /home/local/USHERBROOKE/ouew2201/Documents/GenePredictionReviewBenchmark/Test/Predictions/helixer/150bp/helixer_A4I4C4_LEIIN.fasta.gff3")
            os.system("python3.10 Helixer/Helixer.py --lineage invertebrate --fasta-path /home/local/USHERBROOKE/ouew2201/Documents/GenePredictionReviewBenchmark/g3po_main/References/Fasta_confirmed/150bp/B3M138_DROAN.fasta --gff-output-path /home/local/USHERBROOKE/ouew2201/Documents/GenePredictionReviewBenchmark/Test/Predictions/helixer/150bp/helixer_B3M138_DROAN.fasta.gff3")
            
            return True
        except:
            return False
    else:
        try:
            if add_nuc == 150:
                output = os.path.join(work_path, "Predictions", "helixer", "150bp")
                output_protein = work_path + f"Predictions/Proteines/helixer/150bp/"
                fasta_dir = os.path.join(work_path, "References", "Fasta_confirmed", "150bp")

            else:
                output = os.path.join(work_path, "Predictions", "helixer", f"{str(add_nuc)}Kb")
                output_protein = work_path + f"Predictions/Proteines/helixer/{add_nuc}Kb/"
                fasta_dir = os.path.join(work_path, "References", "Fasta_confirmed", f"{str(add_nuc)}Kb/")
            
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
                        #cmd = f"python3 {os.path.join(path_soft, "Helixer.py")} --lineage {species} --fasta-path {os.path.join(fasta_dir, elmt)} --gff-output-path {os.path.join(output, f"helixer_{elmt}.gff3")}"
                        #cmd = cmd = f"python3 {os.path.join(path_soft, 'Helixer.py')} --lineage {species} --fasta-path {os.path.join(fasta_dir, elmt)} --gff-output-path {os.path.join(output, f'helixer_{elmt}.gff3')}"
                        cmd = f"python3.10 Helixer/Helixer.py --lineage {species} --fasta-path {os.path.join(fasta_dir, elmt)} --gff-output-path {os.path.join(output, f'helixer_{elmt}.gff3')}"

                        print(cmd)
                        os.system(cmd)
                #break
            
            t1 = time.time()
            print(f"Helixer took {t1 - t0} seconds for {add_nuc}Kb")
        
        except:
            return False

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
    #time_sortie.write(f"HMMgene\t{add_nuc}\t{t_final}\n")
    print(f"HMMgene took {t_final} seconds ")
    print(f"{add_nuc}Kb")


if __name__ == '__main__':
    # workpath is the root directory where you download the git repository
    workpath =  os.path.join(os.getcwd(), "g3po_main")
    testpath =  os.path.join(os.getcwd(), "Test")

    # Please change the empty fields by the path of the different programs (eg: "/home/user/Genscanlinux/"")
    path_soft_augustus = "/usr/bin"
    path_soft_helixer = os.path.join(os.getcwd(), "Helixer")



    tt0 = time.perf_counter()
    _nuc = [150, 2, 4, 6, 8, 10]
    #time_sortie = open("./time_prediction.csv", "w")

    for add_nuc in _nuc:
        launch_augustus(workpath, add_nuc, path_soft_augustus, True)
        launch_helixer(workpath, add_nuc, path_soft_helixer, True) 
        #launch_hmmgene(workpath, add_nuc, path_soft_hmmgene)
        #break

    #tt1 = time.perf_counter()

    #final_time = tt1-tt0
    
    #time_sortie.close()
    #print(f"time clock : {final_time}")

    #time_sortie.write(f"Helixer\t{add_nuc}\t{tt1 - tt0}\n")