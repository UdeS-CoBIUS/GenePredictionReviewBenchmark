#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created : 25/08/2019

@author: Scalzitti Nicolas
"""

import os


def define_dna_length(line):
    longueur_sequence = line.split("=")
    longueur_sequence = int(longueur_sequence[1].split(",")[0].strip())
    return longueur_sequence


def augustus_to_gar(work_path, program, add_nuc, rem_file=False):
    """
    Converts prediction file from Augustus to .gar file
    :param work_path: chemin du dossier de travail (root)
    :param program: Name of the predictor
    :param add_nuc: Number of flanked DNA
    :param rem_file: if True, remove files without prediction
    :return:
    """

    if add_nuc != 150:
        output = work_path + f"Predictions/Exon_maps/{program}/{add_nuc}Kb/"
        os.makedirs(output, exist_ok=True)

        path_prediction = work_path + f"Predictions/{program}/{add_nuc}Kb/"
        liste_file = os.listdir(path_prediction)
        add_nuc = int(add_nuc) * 1000

    else:
        add_nuc = int(add_nuc)
        output = work_path + f"Predictions/Exon_maps/{program}/150bp/"
        os.makedirs(output, exist_ok=True)
        path_prediction = work_path + f"Predictions/{program}/150bp/"
        liste_file = os.listdir(path_prediction)

    nbr_file_without_prediction = 0

    for file in liste_file:
        with open(path_prediction + file, "r") as file_1:
            for ligne in file_1:
                if "# (none)" in ligne:
                    nbr_file_without_prediction += 1
                    if rem_file is True:
                        os.remove(path_prediction + file)
                    else:
                        pass

        with open(path_prediction + file, "r") as file_X:

            # open each prediction file
            id_uniprot = file.replace("augustus_", "").replace(".fasta", "")
            liste_gene = []
            dico_ligne = {}
            longueur_sequence = 0

            # Each line of the file is sorted, we retrieve items from the list if they contain info from gene strucuture
            for ligne in file_X:

                # Determines the length of the original sequence
                if "length" in ligne:
                    longueur_sequence = define_dna_length(ligne)

                if "#" not in ligne:
                    # We take the existing number of transcripts and add it to the dictionary
                    if 'gene ' in ligne or "gene\t" in ligne:
                        ligne = ligne.split("\t")  # ou split("\t")
                        ligne = list(filter(None, ligne))

                        if ligne[-1].strip() in dico_ligne.keys():
                            dico_ligne[ligne[-1].strip()].append(ligne)

                        else:
                            dico_ligne[ligne[-1].strip()] = []
                            dico_ligne[ligne[-1].strip()].append(ligne)

                        # Contains the line with the gene number on it, ex:
                        # ['BBS9|G1TZP3_RABIT|ENSO', 'AUGUSTUS', 'gene', '6786', '39597', '1', '+', '.', 'g1\n']
                        liste_gene.append(ligne)

            liste_score = []
            dico_longueur_gene = {}

            for k, v, in dico_ligne.items():
                long = (int(v[0][4])) - (int(v[0][3]))
                dico_longueur_gene[k] = [long]

                strand = v[0][6]
                for val in v:
                    if strand == "+":
                        liste_score.append(int(val[3]))
                        dico_longueur_gene[k].append(int(val[3]))
                        dico_longueur_gene[k].append(int(val[4]))

                    elif strand == "-":
                        liste_score.append(int(val[4]))
                        dico_longueur_gene[k].append(int(val[4]))
                        dico_longueur_gene[k].append(int(val[3]))

            maximum = 0
            best_transcrit = "None"

            for k, v in dico_longueur_gene.items():

                # start and stop prediction - the added nucleotides
                start = int(v[1]) - add_nuc
                stop = int(v[2]) - add_nuc

                # Initialization of the boundaries
                borne = [0, 0]

                # Start and stop of the reference
                stop_ref = longueur_sequence - (add_nuc * 2)
                start_ref = 1

                # Prediction is totally ahead of the gene
                if start < 0 and stop < 0:
                    pass
                # The prediction is totally downstream of the gene
                elif start > stop_ref and stop > stop_ref:
                    pass
                # The beginning of the prediction is out of the gene and the rest is in the gene
                elif start < 0 and stop <= stop_ref:
                    borne = [1, stop]
                    pass
                # The beginning and the end are in the gene
                elif start > 0 and stop <= stop_ref:
                    borne = [start, stop]
                    pass
                # The beginning is in the gene and the end is out
                elif start > 0 and stop > stop_ref:
                    borne = [start, stop_ref]
                    pass
                # The beginning and the end overlap the entire gene
                elif start < 0 and stop > stop_ref:
                    borne = [start_ref, stop_ref]
                    pass
                else:
                    pass

                # length of the predicted sequence that corresponds to the reference gene
                long = abs(borne[0] - borne[1])

                # k is the best prediction
                if long > maximum:
                    maximum = long
                    best_transcrit = k

            print(">>> The gene of interest for : " + str(file).replace("augustus_", "").replace(".fasta", "") +
                  " is : " + str(best_transcrit))

            if best_transcrit == "None":
                print(">>> No good prediction in " + file)
                with open(output + file.replace(".fasta", "") + ".gar", "w") as file_E:
                    file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_Augustus\n")
                    file_E.write("> 0 exon\n\n")
            else:

                len_prot = 0

                with open(path_prediction + file, "r") as file_X2:
                    liste_cds = []
                    liste_exon = []
                    start_codon = "None"
                    stop_codon = "None"
                    tss = "None"
                    tts = "None"

                    for i, ligne in enumerate(file_X2):
                        if "#" in ligne:
                            pass
                        else:
                            if best_transcrit in ligne and "gene\t" in ligne:
                                ligne = ligne.strip().split("\t")
                                strand = ligne[6]

                                if strand == "-":
                                    start2 = int(ligne[4]) - add_nuc
                                    stop2 = int(ligne[3]) - add_nuc

                                elif strand == "+":

                                    start2 = int(ligne[3]) - add_nuc
                                    stop2 = int(ligne[4]) - add_nuc

                                len_prot = abs((stop2 - start2) + 1)


                            if "exon\t" in ligne and best_transcrit in ligne:
                                liste_exon.append(ligne.strip().split("\t"))

                            # TSS
                            if "tss\t" in ligne and best_transcrit in ligne:
                                tss = int(ligne.strip().split("\t")[3]) - add_nuc

                            # TTS
                            if "tts\t" in ligne and best_transcrit in ligne:
                                tts = int(ligne.strip().split("\t")[4]) - add_nuc

                            # We add to a list all the lines that contain CDSs
                            if "CDS\t" in ligne and best_transcrit in ligne:
                                liste_cds.append(ligne.strip().split("\t"))
                            # ATG
                            if "start_codon" in ligne and best_transcrit in ligne:
                                start_codon = int(ligne.strip().split("\t")[3]) - add_nuc

                            #Codon Stop
                            if "stop_codon" in ligne and best_transcrit in ligne:
                                stop_codon = int(ligne.strip().split("\t")[4]) - add_nuc

                    liste_myformat_cds = []
                    liste_myformat_exon = []
                    len_cds = 0
                    strand = liste_cds[0][6]

                    if strand == "+":
                        strand = "1"
                    if strand == "-":
                        strand = "-1"

                    # Exon
                    #if len(liste_exon) >0:
                    maximum = len(liste_exon)

                    for i, ligne in enumerate(liste_exon):
                        long = (int(ligne[4]) - (start - 1)) - (int(ligne[3]) - (start - 1))
                        # 1 exon
                        new_start = (int(ligne[3]) - add_nuc)
                        new_stop = (int(ligne[4]) - add_nuc)

                        if maximum == 1:
                            if tss != "None":
                                if int(ligne[3]) - add_nuc != tss:
                                    new_start = (int(ligne[3]) - add_nuc)

                            if tts != "None":
                                if int(ligne[4]) - add_nuc != tts:
                                    new_stop = (int(ligne[4]) - add_nuc)

                            line = "Exon_" + str(i + 1) + "\t" + str(new_start) + "\t" + \
                                   str(new_stop) + "\t" + str(abs(new_stop - new_start) + 1) + "\t" + str(
                                strand) + "\n"
                        else:
                            # 2 or more exons
                            # exon start
                            if i == 0 and tss != "None":
                                if int(ligne[3]) - add_nuc != tss:
                                    new_start = (int(ligne[3]) - add_nuc) - 3
                                    line = "Exon_" + str(i + 1) + "\t" + str(new_start) + "\t" + \
                                           str((int(ligne[4])) - add_nuc) + "\t" + str(
                                        abs((int(ligne[4])) - add_nuc) - new_start + 1) + "\t" + str(
                                        strand) + "\n"
                                else:
                                    line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[3]) - add_nuc)) + "\t" + \
                                           str((int(ligne[4])) - add_nuc) + "\t" + str(long + 1) + "\t" + str(
                                        strand) + "\n"

                            # exon stop
                            elif i == maximum - 1 and tts != "None":
                                if int(ligne[4]) - add_nuc != tts:
                                    new_stop = (int(ligne[4]) - add_nuc) + 3
                                    line = "Exon_" + str(i + 1) + "\t" + str(int(ligne[3]) - add_nuc) + "\t" + \
                                           str(new_stop) + "\t" + str(
                                        abs(new_stop - (int(ligne[3]) - add_nuc)) + 1) + "\t" + str(
                                        strand) + "\n"
                                else:
                                    line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[3]) - add_nuc)) + "\t" + \
                                           str((int(ligne[4])) - add_nuc) + "\t" + str(long + 1) + "\t" + str(
                                        strand) + "\n"

                            # internal exon
                            else:
                                line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[3]) - add_nuc)) + "\t" + \
                                       str((int(ligne[4])) - add_nuc) + "\t" + str(long + 1) + "\t" + str(
                                    strand) + "\n"

                        # len_cds += long + 1
                        liste_myformat_exon.append(line)

                    # CDS
                    maximum = len(liste_cds)

                    for i, ligne in enumerate(liste_cds):
                        long = (int(ligne[4]) - (start - 1)) - (int(ligne[3]) - (start - 1))
                        # 1 exon
                        new_start = (int(ligne[3]) - add_nuc)
                        new_stop = (int(ligne[4]) - add_nuc)

                        if maximum == 1:
                            if start_codon != "None":
                                if int(ligne[3]) - add_nuc != start_codon:
                                    new_start = (int(ligne[3]) - add_nuc) - 3

                            if stop_codon != "None":
                                if int(ligne[4]) - add_nuc != stop_codon:
                                    new_stop = (int(ligne[4]) - add_nuc) + 3

                            line = "Exon_" + str(i + 1) + "\t" + str(new_start) + "\t" + \
                                   str(new_stop) + "\t" + str(abs(new_stop - new_start) + 1) + "\t" + str(
                                strand) + "\n"
                        else:
                            # 2 or more exons
                            # exon start
                            if i == 0 and start_codon != "None":
                                if int(ligne[3]) - add_nuc != start_codon:
                                    new_start = (int(ligne[3]) - add_nuc) - 3
                                    line = "Exon_" + str(i + 1) + "\t" + str(new_start) + "\t" + \
                                           str((int(ligne[4])) - add_nuc) + "\t" + str(
                                        abs((int(ligne[4])) - add_nuc) - new_start + 1) + "\t" + str(
                                        strand) + "\n"
                                else:
                                    line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[3]) - add_nuc)) + "\t" + \
                                           str((int(ligne[4])) - add_nuc) + "\t" + str(long + 1) + "\t" + str(
                                        strand) + "\n"

                            # exon stop
                            elif i == maximum - 1 and stop_codon != "None":
                                if int(ligne[4]) - add_nuc != stop_codon:
                                    new_stop = (int(ligne[4]) - add_nuc) + 3
                                    line = "Exon_" + str(i + 1) + "\t" + str(int(ligne[3]) - add_nuc) + "\t" + \
                                           str(new_stop) + "\t" + str(
                                        abs(new_stop - (int(ligne[3]) - add_nuc)) + 1) + "\t" + str(
                                        strand) + "\n"
                                else:
                                    line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[3]) - add_nuc)) + "\t" + \
                                           str((int(ligne[4])) - add_nuc) + "\t" + str(long + 1) + "\t" + str(
                                        strand) + "\n"

                            # internal exon
                            else:
                                line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[3]) - add_nuc)) + "\t" + \
                                       str((int(ligne[4])) - add_nuc) + "\t" + str(long + 1) + "\t" + str(
                                    strand) + "\n"

                        liste_myformat_cds.append(line)

                    #File CDS
                    with open(output + file.replace(".fasta", "") + "_cds.gar", "w") as file_E:
                        file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_Augustus(CDS)\n")
                        file_E.write("> " + str(len(liste_cds)) + " exons\n\n")

                        for ligne in liste_myformat_cds:
                            file_E.write(ligne)

                    # File exon
                    if len(liste_myformat_exon) == 0:
                       pass
                    else:
                        with open(output + file.replace(".fasta", "") + "_exon.gar", "w") as file_E:
                            file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_Augustus(exon)\n")
                            file_E.write("> " + str(len(liste_exon)) + " exons\n\n")
                            for ligne in liste_myformat_exon:
                                file_E.write(ligne)

    print(f"\n>>> {str(nbr_file_without_prediction)} file(s) have no predictions")


def genscan_to_gar(work_path, program, add_nuc, rem_file=False):
    """
    Converts prediction file from Genscan to .gar file
    :param work_path: chemin du dossier de travail (root)
    :param program: Name of the predictor
    :param add_nuc: Number of flanked DNA
    :param rem_file: if True, remove files without prediction
    :return:
    """

    if add_nuc != 150:
        output = work_path + f"Predictions/Exon_maps/{program}/{add_nuc}Kb/"
        os.makedirs(output, exist_ok=True)

        path_prediction = work_path + f"Predictions/{program}/{add_nuc}Kb/"
        liste_file = os.listdir(path_prediction)
        add_nuc = int(add_nuc) * 1000

    else:
        add_nuc = int(add_nuc)
        output = work_path + f"Predictions/Exon_maps/{program}/150bp/"
        os.makedirs(output, exist_ok=True)
        path_prediction = work_path + f"Predictions/{program}/150bp/"
        liste_file = os.listdir(path_prediction)

    nbr_file_without_prediction = 0
    no_good_pred = 0

    for file in liste_file:
        with open(path_prediction + file, "r") as file_1:
            for ligne in file_1:
                if "NO EXONS/GENES PREDICTED IN SEQUENCE" in ligne:
                    nbr_file_without_prediction += 1

        # Predictions ok
        with open(path_prediction + file, "r") as file_X:
            id_uniprot = file.replace("genscan_", "").replace(".fasta", "")

            liste_gene = []
            dico_ligne = {}
            longueur_sequence = 0

            for ligne in file_X:

                if "bp" in ligne:
                    longueur_sequence = ligne.split(":")
                    longueur_sequence = int(longueur_sequence[1].split(" ")[1].strip())

                if 'Intr' in ligne or "Term" in ligne or "Sngl" in ligne or \
                        "Init" in ligne:
                    ligne = ligne.split(" ")
                    ligne = list(filter(None, ligne))

                    num_transcrit = ligne[0].split(".")[0]

                    if num_transcrit in dico_ligne.keys():
                        dico_ligne[num_transcrit].append(ligne)
                    else:
                        dico_ligne[num_transcrit] = []
                        dico_ligne[num_transcrit].append(ligne)

                    liste_gene.append(ligne)

            liste_score = []

            dico_longueur_gene = {}
            dico_start_exon = {}
            dico_stop_exon = {}

            for k, v, in dico_ligne.items():
                dico_longueur_gene[k] = []
                dico_start_exon[k] = []
                dico_stop_exon[k] = []

                strand = v[0][2]
                for val in v:
                    if val[1] == 'Intr' or val[1] == "Term" or val[1] == "Sngl" or val[1] == "Init":

                        if strand == "+":
                            liste_score.append(int(val[3]))
                            dico_start_exon[k].append(int(val[3]))
                            dico_stop_exon[k].append(int(val[4]))

                        elif strand == "-":
                            liste_score.append(int(val[4]))
                            dico_start_exon[k].append(int(val[4]))
                            dico_stop_exon[k].append(int(val[3]))

            dico_borne = {}
            for k in dico_start_exon.keys():
                dico_borne[k] = []

            for num_transcrit, borne_start in dico_start_exon.items():
                dico_borne[num_transcrit].append(min(borne_start))

            for num_transcrit, borne_stop in dico_stop_exon.items():
                dico_borne[num_transcrit].append(max(borne_stop))

            maximum = 0
            best_transcrit = "None"

            for k, v in dico_borne.items():

                start = int(v[0]) - add_nuc
                stop = int(v[1]) - add_nuc
                borne = [0, 0]

                stop_ref = longueur_sequence - (add_nuc * 2)
                start_ref = 1

                if start < 0 and stop < 0:
                    pass
                elif start > stop_ref and stop > stop_ref:
                    pass
                elif start < 0 and stop <= stop_ref:
                    borne = [1, stop]
                    pass
                elif start > 0 and stop <= stop_ref:
                    borne = [start, stop]
                    pass
                elif start > 0 and stop > stop_ref:
                    borne = [start, stop_ref]
                    pass
                elif start < 0 and stop > stop_ref:
                    borne = [start_ref, stop_ref]
                    pass
                else:
                    pass

                long2 = abs(borne[0] - borne[1])

                if long2 > maximum:
                    maximum = long2
                    best_transcrit = k

            print(">>> The gene of interest for : " + str(file).replace("genscan_", "").replace(".fasta", "") +
                  " is : " + str(best_transcrit))

            if best_transcrit == "None":
                print(">>> No good prediction in " + file)
                no_good_pred+=1
                with open(output + file.replace(".fasta", "") + ".gar", "w") as file_E:
                    file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_Genscan\n")
                    file_E.write("> 0 exon\n\n")
            else:
                with open(path_prediction + file, "r") as file_X2:
                    liste_cds = []
                    cds = 0

                    for i, ligne in enumerate(file_X2):
                        line = ligne.strip().split(" ")[0].split(".")[0]

                        if line != best_transcrit:
                            pass
                        else:
                            if 'Intr' in ligne or "Term" in ligne or "Sngl" in ligne or "Init" in ligne:
                                liste_cds.append(list(filter(None, ligne.strip().split(" "))))

                                ligne = ligne.strip().split(" ")
                                ligne = list(filter(None, ligne))
                                cds += int(ligne[5])

                    liste_myformat = []
                    strand = liste_cds[0][2]

                    if strand == "+":
                        for i, ligne in enumerate(liste_cds):
                            long = int(ligne[5])
                            strand = "1"
                            line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[3]) - add_nuc)) + "\t" + \
                                   str((int(ligne[4])) - add_nuc) + "\t" + str(long) + "\t" + str(strand) + "\n"
                            liste_myformat.append(line)

                    elif strand == "-":
                        for i, ligne in enumerate(reversed(liste_cds)):
                            long = int(ligne[5])
                            strand = "-1"
                            line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[4]) - add_nuc)) + "\t" + \
                                   str((int(ligne[3])) - add_nuc) + "\t" + str(long) + "\t" + str(strand) + "\n"
                            liste_myformat.append(line)

                    with open(output + file.replace(".fasta", "") + ".gar", "w") as file_E:
                        file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_Genscan\n")
                        file_E.write("> " + str(len(liste_cds)) + " exons\n\n")

                        for ligne in liste_myformat:
                            file_E.write(ligne)

    print(f"\n>>> {str(nbr_file_without_prediction)} file(s) have no predictions")
    print(f"\n>>> {str(no_good_pred)} file(s) have no good predictions")


def geneid_to_gar(work_path, program, add_nuc):
    """
    Converts prediction file from GeneID to .gar file
    :param work_path: chemin du dossier de travail (root)
    :param program: Name of the predictor
    :param add_nuc: Number of flanked DNA
    :return:
    """

    if add_nuc != 150:
        output = work_path + f"Predictions/Exon_maps/{program}/{add_nuc}Kb/"
        os.makedirs(output, exist_ok=True)

        path_prediction = work_path + f"Predictions/{program}/{add_nuc}Kb/"
        liste_file = os.listdir(path_prediction)
        add_nuc = int(add_nuc) * 1000

    else:
        add_nuc = int(add_nuc)
        output = work_path + f"Predictions/Exon_maps/{program}/150bp/"
        os.makedirs(output, exist_ok=True)
        path_prediction = work_path + f"Predictions/{program}/150bp/"
        liste_file = os.listdir(path_prediction)

    nbr_file_without_prediction = 0

    for file in liste_file:
        with open(path_prediction + file, "r") as file_X:

            id_uniprot = file.replace("geneid_", "").replace(".fasta", "")
            liste_gene = []
            dico_ligne = {}

            for ligne in file_X:

                if "Length" in ligne:
                    longueur_sequence = ligne.split(" ")
                    longueur_sequence = int(longueur_sequence[6])

                if "#" not in ligne:
                    if 'Internal' in ligne or "Terminal" in ligne or "First" in ligne or 'Single' in ligne:
                        ligne = ligne.split("\t")  # ou split("\t")
                        ligne = list(filter(None, ligne))
                        cle = ligne[-1].strip().split("|")[-1]

                        if cle in dico_ligne.keys():
                            dico_ligne[cle].append(ligne)
                        else:
                            cle = ligne[-1].strip().split("|")[-1]

                            dico_ligne[cle] = []
                            dico_ligne[cle].append(ligne)

                        liste_gene.append(ligne)

            dico_longueur_gene = {}
            dico_start_exon = {}
            dico_stop_exon = {}

            for k, v in dico_ligne.items():
                dico_longueur_gene[k] = []
                dico_start_exon[k] = []
                dico_stop_exon[k] = []
                liste_score = []

                for val in v:
                    bornes = val[0].split(" ")
                    bornes = list(filter(None, bornes))

                    interne = int(bornes[1])
                    externe = int(bornes[2])

                    liste_score.append([int(interne), int(externe)])

                mini = liste_score[0][0]

                for i, elmt in enumerate(liste_score):

                    mini_tpr = min(liste_score[i])
                    if mini_tpr <= mini:
                        mini = mini_tpr

                dico_start_exon[k].append(mini)

                maxi = liste_score[0][0]
                for i, elmt in enumerate(liste_score):

                    maxi_tpr = max(liste_score[i])
                    if maxi_tpr >= maxi:
                        maxi = maxi_tpr

                dico_stop_exon[k].append(maxi)

            dico_borne = {}
            for k in dico_start_exon.keys():
                dico_borne[k] = []

            for num_transcrit, borne_start in dico_start_exon.items():
                dico_borne[num_transcrit].append(min(borne_start))

            for num_transcrit, borne_stop in dico_stop_exon.items():
                dico_borne[num_transcrit].append(max(borne_stop))

            maximum = 0
            best_transcrit = "None"

            for k, v in dico_borne.items():

                start = int(v[0]) - add_nuc
                stop = int(v[1]) - add_nuc

                if start > stop:
                    start = int(v[1]) - add_nuc
                    stop = int(v[0]) - add_nuc

                borne = [0, 0]

                stop_ref = longueur_sequence - (add_nuc * 2)
                start_ref = 1

                if start < 0 and stop < 0:
                    pass
                elif start > stop_ref and stop > stop_ref:
                    pass
                elif start < 0 and stop <= stop_ref:
                    borne = [1, stop]
                    pass
                elif start > 0 and stop <= stop_ref:
                    borne = [start, stop]
                    pass
                elif start > 0 and stop > stop_ref:
                    borne = [start, stop_ref]
                    pass
                elif start < 0 and stop > stop_ref:
                    borne = [start_ref, stop_ref]
                    pass
                else:
                    pass

                long2 = abs(borne[0] - borne[1]) + 1

                if long2 > maximum:
                    maximum = long2
                    best_transcrit = k

        print(">>> The gene of interest for : " + str(file).replace("geneid_", "").replace(".fasta", "") +
              " is : " + str(best_transcrit))
        if best_transcrit == "None":
            print(">>> No good prediction in " + file)
            with open(output + file.replace(".fasta", "") + ".gar", "w") as file_E:
                file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_GeneID\n")
                file_E.write("> 0 exon\n\n")
        else:
            with open(path_prediction + file, "r") as file_X2:
                liste_cds = []
                cds = 0

                for i, ligne in enumerate(file_X2):
                    if "#" not in ligne:
                        line = ligne.strip().split("\t")[-1].split("|")[-1]

                        if line != best_transcrit:
                            pass
                        else:
                            if 'Internal' in ligne or "Terminal" in ligne or "First" in ligne or 'Single' in ligne:
                                ligne = ligne.strip().split("\t")
                                ligne = list(filter(None, ligne))
                                liste_cds.append(ligne)
                                ligne = ligne[0].split(" ")
                                ligne = list(filter(None, ligne))

                                a = int(ligne[1])
                                b = int(ligne[2])
                                cds += abs(a - b)

                    else:
                        pass

                liste_myformat = []
                strand = liste_cds[0][2].split(" ")[0]

                for i, ligne in enumerate(liste_cds):

                    ligne = ligne[0].split(" ")
                    ligne = list(filter(None, ligne))

                    a = int(ligne[1])
                    b = int(ligne[2])
                    long = abs(a - b)

                    if strand == "+":
                        strand = "1"
                    elif strand == "-":
                        strand = "-1"

                    line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[1]) - add_nuc)) + "\t" + \
                           str((int(ligne[2])) - add_nuc) + "\t" + str(long) + "\t" + str(strand) + "\n"
                    liste_myformat.append(line)

                for i, elmt in enumerate(liste_myformat):
                    if i == 0:
                        a = elmt.strip().split("\t")[1]
                    if i + 1 == len(liste_cds):
                        b = elmt.strip().split("\t")[2]

                len_prot = abs(int(a) - int(b)) + 1

                with open(output + file.replace(".fasta", "") + ".gar", "w") as file_E:
                    file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_GeneID\n")
                    file_E.write("> " + str(len(liste_cds)) + " exons\n\n")

                    for ligne in liste_myformat:
                        file_E.write(ligne)

    print(f"\n>>> {str(nbr_file_without_prediction)} file(s) have no predictions")


def glimmer_to_gar(work_path, program, add_nuc, rem_file=False):
    """
    Converts prediction file from GlimmerHMM to .gar file
    :param work_path: chemin du dossier de travail (root)
    :param program: Name of the predictor
    :param add_nuc: Number of flanked DNA
    :param rem_file: if True, remove files without prediction
    :return:
    """

    if add_nuc != 150:
        output = work_path + f"Predictions/Exon_maps/{program}/{add_nuc}Kb/"
        os.makedirs(output, exist_ok=True)

        path_prediction = work_path + f"Predictions/{program}/{add_nuc}Kb/"
        liste_file = os.listdir(path_prediction)
        add_nuc = int(add_nuc) * 1000

    else:
        add_nuc = int(add_nuc)
        output = work_path + f"Predictions/Exon_maps/{program}/150bp/"
        os.makedirs(output, exist_ok=True)
        path_prediction = work_path + f"Predictions/{program}/150bp/"
        liste_file = os.listdir(path_prediction)

    nbr_file_without_prediction = 0

    for file in liste_file:

        size = os.path.getsize(path_prediction + file)
        if size < 200:
            activator = True
            nbr_file_without_prediction += 1
            if rem_file is True:
                os.remove(path_prediction + file)

        # Predictions ok
        with open(path_prediction + file, "r") as file_X:
            id_uniprot = file.replace("glimmer_", "").replace(".fasta", "")

            liste_gene = []
            dico_ligne = {}
            longueur_sequence = 0

            for ligne in file_X:

                if "##sequence" in ligne:
                    longueur_sequence = ligne.split("|")[-1]
                    longueur_sequence = int(longueur_sequence.strip().split(" ")[2])

                if "CDS\t" in ligne:
                    ligne = ligne.split("\t")
                    ligne = list(filter(None, ligne))

                    num_transcrit = ligne[8].split("|")[-1].split(".")[2].split(";")[0]

                    if num_transcrit in dico_ligne.keys():
                        dico_ligne[num_transcrit].append(ligne)
                    else:
                        dico_ligne[num_transcrit] = []
                        dico_ligne[num_transcrit].append(ligne)

                    liste_gene.append(ligne)

            liste_score = []

            dico_longueur_gene = {}
            dico_start_exon = {}
            dico_stop_exon = {}

            for k, v, in dico_ligne.items():
                dico_longueur_gene[k] = []
                dico_start_exon[k] = []
                dico_stop_exon[k] = []

                for val in v:
                    start = int(val[3])
                    stop = int(val[4])

                    if val[2] == 'CDS':
                        liste_score.append(start)
                        dico_start_exon[k].append(start)
                        dico_stop_exon[k].append(stop)

            dico_borne = {}
            for k in dico_start_exon.keys():
                dico_borne[k] = []

            for num_transcrit, borne_start in dico_start_exon.items():
                dico_borne[num_transcrit].append(min(borne_start))

            for num_transcrit, borne_stop in dico_stop_exon.items():
                dico_borne[num_transcrit].append(max(borne_stop))

            maximum = 0
            best_transcrit = "None"

            for nom_gene, bornes in dico_borne.items():

                start = int(bornes[0]) - add_nuc
                stop = int(bornes[1]) - add_nuc
                borne = [0, 0]

                stop_ref = longueur_sequence - (add_nuc * 2)
                start_ref = 1

                if start < 0 and stop < 0:
                    pass
                elif start > stop_ref and stop > stop_ref:
                    pass
                elif start < 0 and stop <= stop_ref:
                    borne = [1, stop]
                    pass
                elif start > 0 and stop <= stop_ref:
                    borne = [start, stop]
                    pass
                elif start > 0 and stop > stop_ref:
                    borne = [start, stop_ref]
                    pass
                elif start < 0 and stop > stop_ref:
                    borne = [start_ref, stop_ref]
                    pass
                else:
                    pass

                long2 = abs(borne[0] - borne[1])

                if long2 > maximum:
                    maximum = long2
                    best_transcrit = nom_gene

            print(">>> The gene of interest for : " + str(file).replace("genscan_", "").replace(".fasta", "") +
                  " is : " + str(best_transcrit))

            if best_transcrit == "None":
                print(">>> No good prediction in " + file)
                with open(output + file.replace(".fasta", "") + ".gar", "w") as file_E:
                    file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_Glimmer\n")
                    file_E.write("> 0 exon\n\n")

            else:
                with open(path_prediction + file, "r") as file_X2:
                    liste_cds = []
                    cds = 0

                    for i, ligne in enumerate(file_X2):
                        ligne = ligne.split("\t")
                        ligne = list(filter(None, ligne))

                        try:
                            num_tran = ligne[8].split("|")[-1].split(".")[2].split(";")[0]
                        except IndexError:
                            num_tran = ""
                            pass

                        if num_tran != best_transcrit:
                            pass
                        else:
                            if 'CDS' in ligne:
                                liste_cds.append(ligne)

                                start1 = int(ligne[3])
                                stop1 = int(ligne[4])

                                cds += abs(stop1 - start1)
                                strand = ligne[6]

                    cds += 1
                    liste_myformat = []

                    for i, ligne in enumerate(liste_cds):
                        longe = abs(int(ligne[4]) - int(ligne[3])) + 1

                        if strand == "+":
                            strand = "1"
                        elif strand == "-":
                            strand = "-1"
                        line = f"Exon_{str(i + 1)}\t{str((int(ligne[3]) - add_nuc))}\t{str((int(ligne[4])) - add_nuc)}\t{str(longe)}\t{str(strand)}\n"
                        liste_myformat.append(line)

                    with open(output + file.replace(".fasta", "") + ".gar", "w") as file_E:
                        file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_Glimmer\n")
                        file_E.write("> " + str(len(liste_cds)) + " exons\n\n")

                        for ligne in liste_myformat:
                            file_E.write(ligne)

    print(f"\n>>> {str(nbr_file_without_prediction)} file(s) have no predictions")


def snap_to_gar(work_path, program, add_nuc, rem_file=False):
    """
        Converts prediction file from Snap to .gar file
        :param work_path: chemin du dossier de travail (root)
        :param program: Name of the predictor
        :param add_nuc: Number of flanked DNA
        :param rem_file: if True, remove files without prediction
        :return:
        """

    if add_nuc != 150:
        output = work_path + f"Predictions/Exon_maps/{program}/{add_nuc}Kb/"
        os.makedirs(output, exist_ok=True)
        path_prediction = work_path + f"Predictions/{program}/{add_nuc}Kb/"
        liste_file = os.listdir(path_prediction)
        path_fasta = work_path + f"References/Fasta/{add_nuc}Kb/"
        add_nuc = int(add_nuc) * 1000

    else:
        add_nuc = int(add_nuc)
        output = work_path + f"Predictions/Exon_maps/{program}/150bp/"
        os.makedirs(output, exist_ok=True)
        path_prediction = work_path + f"Predictions/{program}/150bp/"
        path_fasta = work_path + "References/Fasta/150bp/"
        liste_file = os.listdir(path_prediction)

    nbr_file_without_prediction = 0

    for file in liste_file:
        size = os.path.getsize(path_prediction + file)
        if size < 10:
            nbr_file_without_prediction += 1
            if rem_file is True:
                os.remove(path_prediction + file)

        # Predictions ok
        with open(path_fasta + file.replace("snap_", ""), "r") as file_R11:
            for ligne in file_R11:
                if ">" not in ligne:
                    longueur_sequence = len(ligne)

        with open(path_prediction + file, "r") as file_X:
            id_uniprot = file.replace("snap_", "").replace(".fasta", "")
            liste_gene = []
            dico_ligne = {}

            for ligne in file_X:

                if 'Exon' in ligne or "Eterm" in ligne or "Einit" in ligne or "Esngl" in ligne:
                    ligne = ligne.strip().split("\t")

                    num_transcrit = ligne[8].split("|")[-1].split(".")[-1]

                    if num_transcrit in dico_ligne.keys():
                        dico_ligne[num_transcrit].append(ligne)
                    else:
                        dico_ligne[num_transcrit] = []
                        dico_ligne[num_transcrit].append(ligne)

                    liste_gene.append(ligne)

            liste_score = []
            dico_longueur_gene = {}
            dico_start_exon = {}
            dico_stop_exon = {}

            for k, v, in dico_ligne.items():
                dico_longueur_gene[k] = []
                dico_start_exon[k] = []
                dico_stop_exon[k] = []

                for val in v:
                    if val[2] == 'Exon' or val[2] == "Eterm" or val[2] == "Einit" or val[2] == "Esngl":
                        liste_score.append(int(val[3]))
                        dico_start_exon[k].append(int(val[3]))
                        dico_stop_exon[k].append(int(val[4]))

            dico_borne = {}
            for k in dico_start_exon.keys():
                dico_borne[k] = []

            for num_transcrit, borne_start in dico_start_exon.items():
                dico_borne[num_transcrit].append(min(borne_start))

            for num_transcrit, borne_stop in dico_stop_exon.items():
                dico_borne[num_transcrit].append(max(borne_stop))

            maximum = 0
            best_transcrit = "None"

            for k, v in dico_borne.items():
                start = int(v[0]) - add_nuc
                stop = int(v[1]) - add_nuc
                borne = [0, 0]
                stop_ref = longueur_sequence - (add_nuc * 2)
                start_ref = 1

                if start < 0 and stop < 0:
                    pass
                elif start > stop_ref and stop > stop_ref:
                    pass
                elif start < 0 and stop <= stop_ref:
                    borne = [1, stop]
                    pass
                elif start > 0 and stop <= stop_ref:
                    borne = [start, stop]
                    pass
                elif start > 0 and stop > stop_ref:
                    borne = [start, stop_ref]
                    pass
                elif start < 0 and stop > stop_ref:
                    borne = [start_ref, stop_ref]
                    pass
                else:
                    pass

                long2 = abs(borne[0] - borne[1])

                if long2 > maximum:
                    maximum = long2
                    best_transcrit = k

            print(">>> The gene of interest for : " + str(file).replace("snap_", "").replace(".fasta", "") +
                  " is : " + str(best_transcrit))

            if best_transcrit == "None":
                print(">>> No good prediction in " + file)
                with open(output + file.replace(".fasta", "") + ".gar", "w") as file_E:
                    file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_Snap\n")
                    file_E.write("> 0 exon\n\n")

            else:
                with open(path_prediction + file, "r") as file_X2:
                    liste_cds = []
                    cds = 0

                    for i, ligne in enumerate(file_X2):
                        num_transcr = ligne.strip().split("\t")[8].split("|")[-1].split(".")[-1]

                        if num_transcr != best_transcrit:
                            pass
                        else:
                            if 'Exon' in ligne or "Eterm" in ligne or "Einit" in ligne or "Esngl" in ligne:
                                ligne = ligne.strip().split("\t")

                                liste_cds.append(list(filter(None, ligne)))
                                cds += abs(int(ligne[4]) - int(ligne[3]))

                    liste_myformat = []
                    strand = liste_cds[0][6]

                    for i, ligne in enumerate(liste_cds):

                        long = abs(int(ligne[4]) - int(ligne[3]))

                        if strand == "+":
                            strand = "1"
                        elif strand == "-":
                            strand = "-1"

                        line = "Exon_" + str(i + 1) + "\t" + str((int(ligne[3]) - add_nuc)) + "\t" + \
                               str((int(ligne[4])) - add_nuc) + "\t" + str(long) + "\t" + str(strand) + "\n"

                        liste_myformat.append(line)

                    liste_start = []
                    liste_stop = []
                    for elmt in liste_myformat:
                        elmt = elmt.strip().split("\t")
                        liste_start.append(elmt[1])
                        liste_stop.append(elmt[2])

                    a = int(min(liste_start))
                    b = int(max(liste_stop))
                    len_mrna = abs(b - a)

                    with open(output + file.replace(".fasta", "") + ".gar", "w") as file_E:
                        file_E.write(f"# Transcrit_{best_transcrit};{id_uniprot};PREDICTED_by_Snap\n")
                        file_E.write("> " + str(len(liste_cds)) + " exons\n\n")

                        for ligne in liste_myformat:
                            file_E.write(ligne)

    print(f"\n>>> {str(nbr_file_without_prediction)} file(s) have no predictions")


if __name__ == '__main__':
    # workpath is the root directory where you download the git repository
    workpath = "my_path/Benchmark_study/"
    _nuc = [150,2,4,6,8,10]

    for add_nuc in _nuc:
        augustus_to_gar(workpath, "augustus", add_nuc)
        genscan_to_gar(workpath , "genscan", add_nuc)
        geneid_to_gar(workpath, "geneid", add_nuc)
        glimmer_to_gar(workpath, "glimmer", add_nuc)
        snap_to_gar(workpath, "snap", add_nuc)

