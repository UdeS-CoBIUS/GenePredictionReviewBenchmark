#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created : 27/08/2019

@author: Scalzitti Nicolas 
"""

from tqdm import tqdm
import os
from pathlib import Path



def moyenne(tableau):
    """
    Calculate the average of the values in array
    :param tableau: an array containing values
    :return:
    """
    return sum(tableau) / len(tableau)


def create_prediction(work_path, soft, add_nuc, file):
    """
    Create the architecture of the prediction file
    :param work_path: The root of the main project
    :param soft: Name of the predictor
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :param file: name of the file .gar
    :return:
    """

    if soft == "augustus":
        dico_transcrit_exon_pred_exon = {}
        dico_transcrit_exon_pred_cds = {}

        if add_nuc == 150:
            file_cds = Path(work_path + f"Predictions/Exon_maps/{soft}/150bp/{file.replace('_exon', '_cds')}")
            file_exon = Path(work_path + f"Predictions/Exon_maps/{soft}/150bp/{file.replace('_cds','_exon')}")
        else:
            file_cds = Path(work_path + f"Predictions/Exon_maps/{soft}/{add_nuc}Kb/{file.replace('_exon', '_cds')}")
            file_exon = Path(work_path + f"Predictions/Exon_maps/{soft}/{add_nuc}Kb/{file.replace('_cds','_exon')}")

        # CDS
        try:
            with open(file_cds, "r") as file_R1:
                for ligne in file_R1:
                    if "Exon_" in ligne:
                        ligne = ligne.strip().split("\t")

                        cle = ligne[0]
                        dico_transcrit_exon_pred_cds[cle] = set()

                        start = int(ligne[1])
                        stop = int(ligne[2])

                        if start >= stop:
                            for nuc in range(stop, start + 1):
                                if nuc == 0:
                                    pass
                                else:
                                    dico_transcrit_exon_pred_cds[cle].add(nuc)

                        elif stop > start:
                            for nuc in range(start, stop + 1):
                                if nuc == 0:
                                    pass
                                else:
                                    dico_transcrit_exon_pred_cds[cle].add(nuc)
        except FileNotFoundError:
            pass

        # EXON
        try:
            with open(file_exon, "r") as file_R2:
                for ligne in file_R2:
                    if "Exon_" in ligne:
                        ligne = ligne.strip().split("\t")

                        cle = ligne[0]
                        dico_transcrit_exon_pred_exon[cle] = set()

                        start = int(ligne[1])
                        stop = int(ligne[2])

                        if start >= stop:
                            for nuc in range(stop, start + 1):
                                if nuc == 0:
                                    pass
                                else:
                                    dico_transcrit_exon_pred_exon[cle].add(nuc)

                        elif stop > start:
                            for nuc in range(start, stop + 1):
                                if nuc == 0:
                                    pass
                                else:
                                    dico_transcrit_exon_pred_exon[cle].add(nuc)

            return dico_transcrit_exon_pred_exon, dico_transcrit_exon_pred_cds
        except FileNotFoundError:
            return dico_transcrit_exon_pred_exon, dico_transcrit_exon_pred_cds

    else:
        dico_transcrit_exon_pred = {}
        if add_nuc == 150:
            file = work_path + f"Predictions/Exon_maps/{soft}/150bp/{file}"
        else:
            file = work_path + f"Predictions/Exon_maps/{soft}/{add_nuc}Kb/{file}"

        with open(file, "r") as file_R1:
            for ligne in file_R1:
                if "Exon_" in ligne:
                    ligne = ligne.strip().split("\t")

                    cle = ligne[0]
                    dico_transcrit_exon_pred[cle] = set()

                    start = int(ligne[1])
                    stop = int(ligne[2])

                    if start >= stop:
                        for nuc in range(stop, start + 1):
                            if nuc == 0:
                                pass
                            else:
                                dico_transcrit_exon_pred[cle].add(nuc)

                    elif stop > start:
                        for nuc in range(start, stop + 1):
                            if nuc == 0:
                                pass
                            else:
                                dico_transcrit_exon_pred[cle].add(nuc)

        return dico_transcrit_exon_pred


def create_reference(work_path, file, soft):
    """
    Create the architecture of the prediction file
    :param work_path: The root of the main project
    :param file: name of the file .gar
    :return:
    """

    dico_transcrit_exon_ref = {}

    with open(work_path + f"References/Exon_map/{file.replace(soft+'_', '').replace('_cds','').replace('_exon','')}" , "r") as file_R1:
        for ligne in file_R1:
            if "Exon_" in ligne:
                ligne = ligne.strip().split("\t")

                cle = ligne[0]
                dico_transcrit_exon_ref[cle] = set()

                start = int(ligne[1])
                stop = int(ligne[2])

                if start >= stop:
                    for nuc in range(stop, start + 1):
                        if nuc == 0:
                            pass
                        else:
                            dico_transcrit_exon_ref[cle].add(nuc)

                elif stop > start:
                    for nuc in range(start, stop + 1):
                        if nuc == 0:
                            pass
                        else:
                            dico_transcrit_exon_ref[cle].add(nuc)

    return dico_transcrit_exon_ref, file


def comparison(work_path, soft, add_nuc, file):
    """
    The function compares the structures of the reference genes with the predicted genes and then calculates and stores
    the different metrics
    :param work_path: The root of the main project
    :param soft: Name of the predictor
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :param file: name of the file (specie name)
    :return:
    """
    if add_nuc == 150:
        output_dir = work_path + f"Resultats/exon/{soft}/150bp/"
    else:
        output_dir = work_path + f"Resultats/exon/{soft}/{str(add_nuc)}Kb/"

    os.makedirs(output_dir + "Details/", exist_ok=True)
    os.makedirs(output_dir + "Calculs/", exist_ok=True)
    os.makedirs(output_dir + "Boundaries/", exist_ok=True)

    # generate the prediction structure
    if soft == "augustus":
        dict_exon_pred_exon = create_prediction(work_path, soft, add_nuc, file)[0]
        dict_exon_pred_cds = create_prediction(work_path, soft, add_nuc, file)[1]

    elif soft == "genscan":
        dict_exon_pred = create_prediction(work_path, soft, add_nuc, file)

    elif soft == "geneid":
        dict_exon_pred = create_prediction(work_path, soft, add_nuc, file)

    elif soft == "glimmer":
        dict_exon_pred = create_prediction(work_path, soft, add_nuc, file)

    elif soft == "snap":
        dict_exon_pred = create_prediction(work_path, soft, add_nuc, file)

    else:
        print(">>> Error in program name, please choose : augustus, genscan, geneid, glimmer or snap")
        exit()

    # generate the reference structure
    dict_exon_ref, file = create_reference(work_path, file, soft)

    # returns the higher nucleotide
    all_nuc = set()
    for exon, nuc in dict_exon_ref.items():
        for nucleotide in nuc:
            all_nuc.add(nucleotide)
    out_sup = max(all_nuc)
    out_inf = min(all_nuc)

    # Initialisation of header
    with open(output_dir + f"Details/{file.replace('.gar', '.csv')}", "a") as file_details:

        file_details.write("Correct Ref;Correct Pred;MissingExon;WrongExon;OutExon;Good_Donnor;Good_Acceptor\n")

        with open(output_dir + f"Calculs/{file.replace('.gar', '.csv')}", "a") as file_calcul:
            file_calcul.write("Correct Ref;Correct Pred;MissingExon;WrongExon;OutExon;Sensibility;Specificity;MEsore;"
                              "WEscore;Good_donnor;Good_acceptor;%_good_donnor;%_good_acceptor\n")

            with open(output_dir + f"Boundaries/{file.replace('.gar', '.csv')}", "a") as file_limite:
                file_limite.write("5' boundary;3' boundary\n")

                if soft == "augustus":
                    correct_exon_ref_exon = []  # good exons in reference
                    correct_exon_pred_exon = []  # good exons in prediction
                    wrong_exon_exon = []  # in prediction
                    out_exon_exon = []  # Exon with lower boundaries 1
                    missing_exon_exon = []  # in reference
                    correct_donneur_exon = []
                    correct_accepteur_exon = []

                    correct_exon_ref_cds = []  # good exons in reference
                    correct_exon_pred_cds = []  # good exons in prediction
                    wrong_exon_cds = []  # in prediction
                    out_exon_cds = []  # Exon with lower boundaries 1
                    missing_exon_cds = []  # in reference
                    correct_donneur_cds = []
                    correct_accepteur_cds = []
                    correct_limite_start_exon = []
                    correct_limite_stop_exon = []
                    correct_limite_start_cds = []
                    correct_limite_stop_cds = []

                    # EXON
                    for exon_ref, nuc_ref in dict_exon_ref.items():
                        for exon_pred, nuc_pred in dict_exon_pred_exon.items():

                            start_exon_pred_exon = min(nuc_pred)
                            start_exon_ref_exon = min(nuc_ref)
                            stop_exon_pred_exon = max(nuc_pred)
                            stop_exon_ref_exon = max(nuc_ref)

                            # EXONS CORRECTS
                            if start_exon_ref_exon == start_exon_pred_exon and stop_exon_ref_exon == stop_exon_pred_exon:
                                correct_exon_ref_exon.append(exon_ref)
                                correct_exon_pred_exon.append(exon_pred)

                            # Exons out inferieur
                            elif start_exon_pred_exon < 1 and stop_exon_pred_exon < 1 and exon_pred not in out_exon_exon:
                                out_exon_exon.append(exon_pred)

                            # Exons out superieur
                            elif start_exon_pred_exon > out_sup and stop_exon_pred_exon > out_sup and exon_pred not in out_exon_exon:
                                out_exon_exon.append(exon_pred)

                            # SS Acceptor
                            if start_exon_pred_exon == start_exon_ref_exon:
                                # Not an exon/intron boundaries but the first Nuc = Start
                                if start_exon_pred_exon != out_inf:
                                    correct_accepteur_exon.append([exon_pred, start_exon_pred_exon])
                                else:
                                    correct_limite_start_exon.append([exon_pred, start_exon_pred_exon])


                            # SS Donnor
                            if stop_exon_pred_exon == stop_exon_ref_exon:
                                # Not an exon/intron boundaries but the last Nuc = Stop
                                if stop_exon_pred_exon != out_sup:
                                    correct_donneur_exon.append([exon_pred, stop_exon_pred_exon])
                                else:
                                    correct_limite_stop_exon.append([exon_pred, stop_exon_pred_exon])

                    for exon_ref, nuc_ref in dict_exon_ref.items():
                        for exon_pred, nuc_pred in dict_exon_pred_exon.items():

                            intersection = nuc_ref.intersection(nuc_pred)  # Nucléotide commun au 2 exons (ref et pred)

                            if len(intersection) != 0 and exon_pred not in correct_exon_pred_exon and exon_pred not in wrong_exon_exon:
                                wrong_exon_exon.append(exon_pred)

                            # Exons Wrong
                            if exon_pred not in correct_exon_pred_exon and exon_pred not in out_exon_exon and exon_pred not in wrong_exon_exon:
                                wrong_exon_exon.append(exon_pred)

                            # Exon missing (or superposed)
                            if exon_ref not in correct_exon_ref_exon and exon_ref not in missing_exon_exon:
                                missing_exon_exon.append(exon_ref)

                    # calculates of metrics
                    try:
                        sensi_exon = len(correct_exon_ref_exon) / (len(correct_exon_ref_exon) + len(missing_exon_exon))
                    except ZeroDivisionError:
                        sensi_exon = 0

                    try:
                        speci_exon = len(correct_exon_ref_exon) / (len(correct_exon_pred_exon) + len(wrong_exon_exon))
                    except ZeroDivisionError:
                        speci_exon = 0

                    try:
                        missex_exon = len(missing_exon_exon) / (len(correct_exon_ref_exon) + len(missing_exon_exon))
                    except ZeroDivisionError:
                        missex_exon = 0

                    try:
                        wrongex_exon = len(wrong_exon_exon) / (len(correct_exon_pred_exon) + len(wrong_exon_exon))
                    except ZeroDivisionError:
                        wrongex_exon = 0

                    # calculates the number of total exons taken into account in the prediction (excluding out)
                    nbr_exons_exon = int(len(correct_exon_pred_exon) + len(wrong_exon_exon))
                    try:
                        if nbr_exons_exon-1 == 0:
                            good_don_pourcent_exon = ""
                            good_acc_pourcent_exon = ""
                        elif nbr_exons_exon == -1:
                            good_don_pourcent_exon = 0.0
                            good_acc_pourcent_exon = 0.0
                        else:
                            good_don_pourcent_exon = round((len(correct_donneur_exon) * 100) / (nbr_exons_exon-1), 2)
                            good_acc_pourcent_exon = round((len(correct_accepteur_exon)* 100) / (nbr_exons_exon-1), 2)

                    except ZeroDivisionError:
                        good_don_pourcent_exon = 0.0
                        good_acc_pourcent_exon = 0.0

                    # CDS
                    for exon_ref, nuc_ref in dict_exon_ref.items():
                        for exon_pred, nuc_pred in dict_exon_pred_cds.items():

                            start_exon_pred = min(nuc_pred)
                            start_exon_ref = min(nuc_ref)
                            stop_exon_pred = max(nuc_pred)
                            stop_exon_ref = max(nuc_ref)

                            # EXONS CORRECTS
                            if start_exon_ref == start_exon_pred and stop_exon_ref == stop_exon_pred:
                                correct_exon_ref_cds.append(exon_ref)
                                correct_exon_pred_cds.append(exon_pred)

                            # Exons out inferieur
                            elif start_exon_pred < 1 and stop_exon_pred < 1 and exon_pred not in out_exon_cds:
                                out_exon_cds.append(exon_pred)

                            # Exons out superieur
                            elif start_exon_pred > out_sup and stop_exon_pred > out_sup and exon_pred not in out_exon_cds:
                                out_exon_cds.append(exon_pred)

                            # SS Acceptor
                            if start_exon_pred == start_exon_ref:
                                if start_exon_pred != out_inf:
                                    correct_accepteur_cds.append([exon_pred, start_exon_pred])
                                else:
                                    correct_limite_start_cds.append([exon_pred, start_exon_pred])

                            # SS Donnor
                            if stop_exon_pred == stop_exon_ref:
                                # Not an exon/intron boundaries but the last Nuc = Stop
                                if stop_exon_pred != out_sup:
                                    correct_donneur_cds.append([exon_pred, stop_exon_pred])
                                else:
                                    correct_limite_stop_cds.append([exon_pred, stop_exon_pred])


                    for exon_ref, nuc_ref in dict_exon_ref.items():
                        for exon_pred, nuc_pred in dict_exon_pred_cds.items():

                            intersection = nuc_ref.intersection(nuc_pred)  # Nucléotide commun au 2 exons (ref et pred)

                            # Predicted exons who overlaped the ref exons are false
                            if len(intersection) != 0 and exon_pred not in correct_exon_pred_cds and exon_pred not in wrong_exon_cds:
                                wrong_exon_cds.append(exon_pred)

                            # Exons Wrong
                            if exon_pred not in correct_exon_pred_cds and exon_pred not in out_exon_cds and exon_pred not in wrong_exon_cds:
                                wrong_exon_cds.append(exon_pred)

                            # Exon missing (or superposed)
                            if exon_ref not in correct_exon_ref_cds and exon_ref not in missing_exon_cds:
                                missing_exon_cds.append(exon_ref)

                    # calculates of metrics
                    try:
                        sensi_cds = len(correct_exon_ref_cds) / (len(correct_exon_ref_cds) + len(missing_exon_cds))
                    except ZeroDivisionError:
                        sensi_cds = 0

                    try:
                        speci_cds = len(correct_exon_ref_cds) / (len(correct_exon_pred_cds) + len(wrong_exon_cds))
                    except ZeroDivisionError:
                        speci_cds = 0

                    try:
                        missex_cds = len(missing_exon_cds) / (len(correct_exon_ref_cds) + len(missing_exon_cds))
                    except ZeroDivisionError:
                        missex_cds= 0

                    try:
                        wrongex_cds = len(wrong_exon_cds) / (len(correct_exon_pred_cds) + len(wrong_exon_cds))
                    except ZeroDivisionError:
                        wrongex_cds = 0

                    # calculates the number of total exons taken into account in the prediction (excluding out)
                    nbr_exons_cds = int(len(correct_exon_pred_cds) + len(wrong_exon_cds))
                    try:
                        if nbr_exons_cds-1 == 0:
                            good_don_pourcent_cds = ""
                            good_acc_pourcent_cds = ""
                        elif nbr_exons_cds == -1:
                            good_don_pourcent_cds = 0.0
                            good_acc_pourcent_cds = 0.0
                        else:
                            good_don_pourcent_cds = round((len(correct_donneur_cds) * 100) / (nbr_exons_cds-1), 2)
                            good_acc_pourcent_cds = round((len(correct_accepteur_cds)* 100) / (nbr_exons_cds-1), 2)
                    except ZeroDivisionError:
                        good_don_pourcent_cds = 0.0
                        good_acc_pourcent_cds = 0.0

                    score_exon = (sensi_exon + speci_exon) / 2
                    score_cds = (sensi_cds + speci_cds) / 2

                    if score_cds >= score_exon:
                        # Writing file
                        file_details.write(f"{correct_exon_ref_cds};{correct_exon_pred_cds};{missing_exon_cds};{wrong_exon_cds};{out_exon_cds};"
                                           f"{correct_donneur_cds};{correct_accepteur_cds};CDS\n")

                        file_calcul.write(
                            f"{len(correct_exon_ref_cds)};{len(correct_exon_pred_cds)};{len(missing_exon_cds)};{len(wrong_exon_cds)};"
                            f"{len(out_exon_cds)};{sensi_cds};{speci_cds};{missex_cds};{wrongex_cds};{len(correct_donneur_cds)};"
                            f"{len(correct_accepteur_cds)};{good_don_pourcent_cds};{good_acc_pourcent_cds};CDS\n")

                        file_limite.write(f"{correct_limite_start_cds};{correct_limite_stop_cds}\n")


                    else:
                        # Writing file
                        file_details.write(
                            f"{correct_exon_ref_exon};{correct_exon_pred_exon};{missing_exon_exon};{wrong_exon_exon};{out_exon_exon};"
                            f"{correct_donneur_exon};{correct_accepteur_exon};EXON\n")

                        file_calcul.write(
                            f"{len(correct_exon_ref_exon)};{len(correct_exon_pred_exon)};{len(missing_exon_exon)};{len(wrong_exon_exon)};"
                            f"{len(out_exon_exon)};{sensi_exon};{speci_exon};{missex_exon};{wrongex_exon};{len(correct_donneur_exon)};"
                            f"{len(correct_accepteur_exon)};{good_don_pourcent_exon};{good_acc_pourcent_exon};EXON\n")

                        file_limite.write(f"{correct_limite_start_exon};{correct_limite_stop_exon}\n")


                else:
                    correct_exon_ref = []  # good exons in reference
                    correct_exon_pred = []  # good exons in prediction
                    wrong_exon = []  # in prediction
                    out_exon = []  # Exon with lower boundaries 1
                    missing_exon = []  # in reference
                    correct_donneur = []
                    correct_accepteur = []
                    correct_limite_start = []
                    correct_limite_stop = []

                    for exon_ref, nuc_ref in dict_exon_ref.items():
                        for exon_pred, nuc_pred in dict_exon_pred.items():
                            # print(exon_pred, nuc_pred)

                            start_exon_pred = min(nuc_pred)
                            start_exon_ref = min(nuc_ref)

                            stop_exon_pred = max(nuc_pred)
                            stop_exon_ref = max(nuc_ref)

                            # EXONS CORRECTS
                            if start_exon_ref == start_exon_pred and stop_exon_ref == stop_exon_pred:
                                correct_exon_ref.append(exon_ref)
                                correct_exon_pred.append(exon_pred)

                            # Exons out inferieur
                            elif start_exon_pred < 1 and stop_exon_pred < 1 and exon_pred not in out_exon:
                                out_exon.append(exon_pred)

                            # Exons out superieur
                            elif start_exon_pred > out_sup and stop_exon_pred > out_sup and exon_pred not in out_exon:
                                out_exon.append(exon_pred)

                            # SS Acceptor
                            if start_exon_pred == start_exon_ref:
                                # if exon boundary isn't the first boundar != 1st exon
                                if start_exon_pred != out_inf:
                                    correct_accepteur.append([exon_pred, start_exon_pred])
                                else:
                                    correct_limite_start.append([exon_pred, start_exon_pred])
                                    # pass

                            # SS Donnor - 3' exon ou 5' intron
                            if stop_exon_pred == stop_exon_ref:
                                # if exon boundary isn't the last boundary != last exon
                                if stop_exon_pred != out_sup:
                                    correct_donneur.append([exon_pred, stop_exon_pred])
                                else:
                                    correct_limite_stop.append([exon_pred, stop_exon_pred])
                                    # pass

                    for exon_ref, nuc_ref in dict_exon_ref.items():
                        for exon_pred, nuc_pred in dict_exon_pred.items():

                            # overlaped exons are considered as false
                            intersection = nuc_ref.intersection(nuc_pred)

                            # Predicted exons who overlaped the ref exons are false
                            if len(intersection) != 0 and exon_pred not in correct_exon_pred and exon_pred not in wrong_exon:
                                wrong_exon.append(exon_pred)

                            # Exons Wrong
                            if exon_pred not in correct_exon_pred and exon_pred not in out_exon and exon_pred not in wrong_exon:
                                wrong_exon.append(exon_pred)

                            # Exon missing (or superposed)
                            if exon_ref not in correct_exon_ref and exon_ref not in missing_exon:
                                missing_exon.append(exon_ref)

                    # calculates of metrics
                    try:
                        sensi = len(correct_exon_ref) / (len(correct_exon_ref) + len(missing_exon))
                    except ZeroDivisionError:
                        sensi = 0

                    try:
                        speci = len(correct_exon_ref) / (len(correct_exon_pred) + len(wrong_exon))
                    except ZeroDivisionError:
                        speci = 0

                    try:
                        missex = len(missing_exon) / (len(correct_exon_ref) + len(missing_exon))
                    except ZeroDivisionError:
                        missex = 0

                    try:
                        wrongex = len(wrong_exon) / (len(correct_exon_pred) + len(wrong_exon))
                    except ZeroDivisionError:
                        wrongex = 0

                    # calculates the number of total exons taken into account in the prediction (excluding out)
                    nbr_exons = int(len(correct_exon_pred) + len(wrong_exon))
                    try:
                        if nbr_exons-1 == 0:
                            good_don_pourcent = ""
                            good_acc_pourcent = ""
                        elif nbr_exons == -1:
                            good_don_pourcent = 0.0
                            good_acc_pourcent = 0.0

                        else:
                            good_don_pourcent = round((len(correct_donneur) * 100) / (nbr_exons-1), 2)
                            good_acc_pourcent = round((len(correct_accepteur) * 100) / (nbr_exons-1), 2)
                    except ZeroDivisionError:
                        good_don_pourcent = 0.0
                        good_acc_pourcent = 0.0

                    # Writing file
                    file_details.write(f"{correct_exon_ref};{correct_exon_pred};{missing_exon};{wrong_exon};{out_exon};"
                                       f"{correct_donneur};{correct_accepteur}\n")

                    file_calcul.write(f"{len(correct_exon_ref)};{len(correct_exon_pred)};{len(missing_exon)};{len(wrong_exon)};"
                                      f"{len(out_exon)};{sensi};{speci};{missex};{wrongex};{len(correct_donneur)};"
                                      f"{len(correct_accepteur)};{good_don_pourcent};{good_acc_pourcent}\n")

                    file_limite.write(f"{correct_limite_start};{correct_limite_stop}\n")

    if soft == "augustus":
        try:
            if score_cds >= score_exon:
                os.remove(output_dir + f"Details/{file.replace('_cds','').replace('_exon','').replace('.gar', '')}_exon.csv")
                os.remove(output_dir + f"Calculs/{file.replace('_cds','').replace('_exon','').replace('.gar', '')}_exon.csv")
                os.remove(output_dir + f"Boundaries/{file.replace('_cds','').replace('_exon','').replace('.gar', '')}_exon.csv")

            else:
                os.remove(output_dir + f"Details/{file.replace('_cds','').replace('_exon','').replace('.gar', '')}_cds.csv")
                os.remove(output_dir + f"Calculs/{file.replace('_cds','').replace('_exon','').replace('.gar', '')}_cds.csv")
                os.remove(output_dir + f"Boundaries/{file.replace('_cds','').replace('_exon','').replace('.gar', '')}_cds.csv")

        except FileNotFoundError:
            pass


def average(work_path, soft, add_nuc, file):
    """
    Extract the value for each metric
    :param work_path: the root of the main project
    :param soft: Name of the predictor
    :param add_nuc: number of nucleotides flanked upstream and downstream of gene
    :param file: name of the file .csv in /Calculs directory
    :return:
    """
    if add_nuc == 150:
        path = work_path + f"Resultats/exon/{soft}/150bp/Calculs/{file.replace('.gar', '.csv')}"
    else:
        path = work_path + f"Resultats/exon/{soft}/{add_nuc}Kb/Calculs/{file.replace('.gar', '.csv')}"

    with open(path, "r") as file_R1:
        for i, ligne in enumerate(file_R1):
            if i != 0:
                ligne = ligne.strip().split(";")
                sensi = float(ligne[5])
                speci = float(ligne[6])
                me = float(ligne[7])
                we = float(ligne[8])
                if ligne[11] == "":
                    pc_don = ""
                else:
                    pc_don = float(ligne[11])
                if ligne[12] == "":
                    pc_acc = ""
                else:
                    pc_acc = float(ligne[12])
            else:
                pass

    return sensi, speci, me, we, pc_don, pc_acc


def calcul_average(work_path, soft, add_nuc, liste_gene):
    """
    Calculates the average of all metrics
    :param work_path: the root of the main project
    :param soft: Name of the predictor
    :param add_nuc: number of nucleotides flanked upstream and downstream of gene
    :param liste_gene: a list containing the name of each file
    :return:
    """
    if add_nuc == 150:
        output_dir = work_path + f"Resultats/exon/{soft}/150bp/"
    else:
        output_dir = work_path + f"Resultats/exon/{soft}/{str(add_nuc)}Kb/"
    avg_sensi = []
    avg_speci = []
    avg_me = []
    avg_we = []
    avg_pc_don = []
    avg_pc_acc = []

    a = 0

    for file in liste_gene:
        try:
            avg_sensi.append(average(work_path, soft, add_nuc, file)[0])
            avg_speci.append(average(work_path, soft, add_nuc, file)[1])
            avg_me.append(average(work_path, soft, add_nuc, file)[2])
            avg_we.append(average(work_path, soft, add_nuc, file)[3])

            if average(work_path, soft, add_nuc, file)[4] == "":
                pass
            else:
                avg_pc_don.append(average(work_path, soft, add_nuc, file)[4])

                if average(work_path, soft, add_nuc, file)[5] == "":
                    pass
                else:
                    avg_pc_acc.append(average(work_path, soft, add_nuc, file)[5])
        except FileNotFoundError:
            pass

    m_sensi = round(moyenne(avg_sensi), 3)
    m_speci = round(moyenne(avg_speci), 3)
    m_me = round(moyenne(avg_me), 3)
    m_we = round(moyenne(avg_we), 3)
    m_don = round(moyenne(avg_pc_don), 3)
    m_acc = round(moyenne(avg_pc_acc), 3)

    with open(output_dir + f"results.csv", "w") as file_W2:
        file_W2.write("Sensibility;Specificity;MissingExon;WrongExon;Donnor;Acceptor\n")
        file_W2.write(f"{m_sensi};{m_speci};{m_me};{m_we};{m_don};{m_acc}\n")

def count_number_good_start_stop(workpath, soft, add_nuc, seq_type):
    if add_nuc == 150:
        path_b = f"Resultats/exon/{soft}/150bp/Boundaries/"
    else:
        path_b = f"Resultats/exon/{soft}/{add_nuc}Kb/Boundaries/"

    liste_validated = []
    with open(f"/home/scalzitti/These/G3PO_Benchmark/Benchmark_study/Sequence_types/{seq_type}.seq", "r") as file_R0:
        for ligne in file_R0:
            ligne = ligne.strip()
            liste_validated.append(ligne)

    nbr_5prime = 0
    nbr_3prime = 0

    for file in liste_validated:
        if soft == "augustus":
            try:
                path = Path(f"{workpath}{path_b}{soft}_{file}_exon.csv")
                open(path)
            except FileNotFoundError:
                try:
                    path = Path(f"{workpath}{path_b}{soft}_{file}_cds.csv")
                    open(path)
                except FileNotFoundError:
                    path = Path(f"{workpath}{path_b}{soft}_{file}.csv")

        else:
            path = Path(f"{workpath}{path_b}{soft}_{file}.csv")

        with open(f"{path}", "r") as file_R1:
            for i, ligne in enumerate(file_R1):
                if i != 0:
                    ligne = ligne.strip().split(";")
                    _5prime = ligne[0] # Start
                    _3prime = ligne[1] # Stop
                    if _5prime != "[]":
                        nbr_5prime += 1
                    if _3prime != "[]":
                        nbr_3prime += 1

    print(f"\n{soft} have {nbr_5prime} 5' boundaries OK (start),  {round((nbr_5prime*100)/len(liste_validated),2)} %")
    print(f"{soft} have {nbr_3prime} 3' boundaries OK (stop),  {round((nbr_3prime*100)/len(liste_validated),2)} %")


def main(workpath):
    seq_type = "Confirmed"
    _nuc = [150, 2, 4, 6, 8, 10]
    l_soft = ["augustus", "genscan", "geneid", "glimmer", "snap"]

    for soft in l_soft:
        for add_nuc in _nuc:

            if add_nuc == 150:
                ref_dir = workpath + f"Predictions/Exon_maps/{soft}/150bp/"
            else:
                ref_dir = workpath + f"Predictions/Exon_maps/{soft}/{add_nuc}Kb/"

            liste_gene = os.listdir(ref_dir)

            for file in tqdm(liste_gene):
                comparison(workpath, soft, add_nuc, file)

            calcul_average(workpath, soft, add_nuc, liste_gene)
            count_number_good_start_stop(workpath, soft, add_nuc, seq_type)


if __name__ == '__main__':
    # workpath is the root directory where you download the git repository
    workpath = "my_path/Benchmark_study/"
    main(workpath)


