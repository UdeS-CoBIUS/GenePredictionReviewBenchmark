#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created : 25/08/2019

@author: Scalzitti Nicolas 
"""

import os
from tqdm import tqdm
from pathlib import Path


def liste_gene(workpath, soft, add_nuc):
    if add_nuc != 150:
        liste_file = os.listdir(workpath + f"References/Fasta/150bp/")
    else:
        liste_file = os.listdir(workpath + f"References/Fasta/150bp/")

    return liste_file


def sensibility(tp, fn):
    """
    According to Guigo et al., 1996,
    :param tp:  true positives
    :param fn: false negatifs
    :return:
    """
    try:
        return round(tp / (tp + fn), 3)
    except ZeroDivisionError:
        return 0.0


def specificity(tp, fp):
    """
    According to Guigo et al., 1996,
    :param tp:  true positives
    :param fp: false positives
    :return:
    """
    try:
        return round(tp / (tp + fp), 3)
    except ZeroDivisionError:
        return 0.0


def f1_score(tp, fp, fn):
    """
    Calculate the F1 score
    :param tp: true positives
    :param fp: false positives
    :param fn: false negatives
    :return:
    """
    try:
        sensi = sensibility(tp, fn)
        speci = specificity(tp, fp)

        return round(2 * ((sensi * speci) / (sensi + speci)), 3)
    except ZeroDivisionError:
        return 0.0


def moyenne(tableau):
    """
    Calculate the average of the values in array
    :param tableau: an array containing values
    :return:
    """
    return sum(tableau) / len(tableau)


def create_reference(work_path, gene, soft):
    """
    The script generate a dictionnary with the transcrit number as key and 2 lists as values. The list number 1 is all
    exonic nucleotids and the list 2 is all intronic nucleotids.
    :param work_path: working directory path (root)
    :param gene: name of the reference file (XXX.gar in ./Structure_map/References/)
    :param soft name of the soft
    :return:
    """
    path_reference = work_path + f"References/Exon_map/"
    gene = gene.replace(soft + "_","").replace("_cds","").replace("_exon","").replace(".fasta", ".gar")

    file = Path(path_reference + gene)
    if file.exists():

        with open(path_reference + gene, 'r') as file_R1:
            dico_transcrit = {}
            for i, ligne in enumerate(file_R1):

                # Transcrit
                if "#" in ligne:
                    cle = ligne.split(";")[0].replace("# ", "")
                    dico_transcrit[cle] = []
                    nuc_exon_pred = set()
                    activator = True
                    activator_un_exon = False

                # just 1 exon
                if "> 1 exons" in ligne:
                    activator_un_exon = True

                # mutliple exons
                if "Exon" in ligne:
                    line = ligne.strip().split("\t")
                    start = int(line[1])
                    stop = int(line[2])

                    if start >= stop:
                        for nucleotide in range(stop, start + 1):
                            nuc_exon_pred.add(nucleotide)
                    elif start < stop:
                        for nucleotide in range(start, stop + 1):
                            nuc_exon_pred.add(nucleotide)

                # found the exon
                if activator_un_exon:
                    dico_transcrit[cle].append(nuc_exon_pred)
                    activator_un_exon = False

                # found all exons
                if "Intron" in ligne:
                    if activator:
                        dico_transcrit[cle].append(nuc_exon_pred)
                        activator = False
                    else:
                        pass

            # add intronic nuc
            for transcrit, values in dico_transcrit.items():
                try:
                    nuc_intron_pred = set()
                    nuc_exonique = values[0]
                    end = max(nuc_exonique)

                    for nuc_intronique in range(1, end + 1):
                        if nuc_intronique not in nuc_exonique:
                            nuc_intron_pred.add(nuc_intronique)

                    dico_transcrit[transcrit].append(nuc_intron_pred)

                except IndexError:
                    # No introns in ref
                    nuc_intron_pred = set()
                    dico_transcrit[transcrit].append(nuc_intron_pred)
                    pass

        return dico_transcrit
    else:
        return {}


def create_prediction_augustus(gene, add_nuc, path_pred):
    """
    Generate the structure of the prediction file from Augustus
    :param gene: name of the file .gar
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :param path_pred: path of the prediction files ./Predictions/{soft}/{value}Kb/
    :return:
    """

    if add_nuc != 150:
        add_nuc = int(add_nuc) * 1000
    else:
        add_nuc = 150

    nuc_intron_pred_exon = set()
    nuc_intron_pred_cds = set()
    nuc_exon_pred_cds = set()
    nuc_exon_pred_exon = set()

    soft = "augustus"
    name = gene.replace(".gar", ".fasta")
    gene = f"{name}"

    file = Path(path_pred + soft + "_" + gene)


    if file.exists():
        a = open(path_pred + soft + "_" + gene, "r")
        # No prediction
        if "# (none)" in a.read():

            return nuc_exon_pred_exon, nuc_exon_pred_cds, nuc_intron_pred_exon, nuc_intron_pred_cds
        else:
            with open(path_pred + soft + "_" + gene, "r") as file_R1:
                # exonic nuc
                for ligne in file_R1:
                    if "CDS\t" in ligne or "start_codon\t" in ligne or "stop_codon\t" in ligne:
                        line = ligne.strip().split("\t")
                        start = int(line[3]) - add_nuc
                        stop = int(line[4]) - add_nuc

                        if start >= stop:
                            for nucleotide in range(stop, start + 1):
                                nuc_exon_pred_cds.add(nucleotide)
                        elif start < stop:
                            for nucleotide in range(start, stop + 1):
                                nuc_exon_pred_cds.add(nucleotide)

                    if "exon\t" in ligne or "tss\t" in ligne or "tts\t" in ligne:
                        line = ligne.strip().split("\t")
                        start = int(line[3]) - add_nuc
                        stop = int(line[4]) - add_nuc

                        if start >= stop:
                            for nucleotide in range(stop, start + 1):
                                nuc_exon_pred_exon.add(nucleotide)
                        elif start < stop:
                            for nucleotide in range(start, stop + 1):
                                nuc_exon_pred_exon.add(nucleotide)

            try:
                # intronic nuc
                # Exon
                try:
                    for nuc_intronique in range(-add_nuc, max(nuc_exon_pred_exon) + add_nuc + 1):
                        if nuc_intronique == 0:
                            pass
                        else:
                            if nuc_intronique not in nuc_intron_pred_exon and nuc_intronique not in nuc_exon_pred_exon:
                                nuc_intron_pred_exon.add(nuc_intronique)

                except ValueError:
                    pass

                # CDS
                try:
                    for nuc_intronique in range(-add_nuc, max(nuc_exon_pred_cds) + add_nuc+ 1):
                        if nuc_intronique == 0:
                            pass
                        else:
                            if nuc_intronique not in nuc_intron_pred_cds and nuc_intronique not in nuc_exon_pred_cds:
                                nuc_intron_pred_cds.add(nuc_intronique)
                except ValueError:
                    pass


            except UnboundLocalError:
                pass

        return nuc_exon_pred_exon, nuc_exon_pred_cds, nuc_intron_pred_exon, nuc_intron_pred_cds
    else:
        print(f"{file} don't exist")
        pass


def create_prediction_genscan(gene, add_nuc, path_pred):
    """
    Generate the structure of the prediction file from Genscan
    :param gene: name of the file .gar
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :param path_pred: path of the prediction files ./Predictions/{soft}/{value}Kb/
    :return:
    """

    soft = "genscan"
    name = gene.replace(".gar", ".fasta")
    gene = f"{name}"

    if add_nuc != 150:
        add_nuc = int(add_nuc) * 1000
    else:
        add_nuc = 150

    file = Path(path_pred + soft + "_" + gene)

    nuc_exon_pred = set()
    nuc_intron_pred = set()

    if file.exists():

        with open(path_pred + soft + "_" + gene, "r") as file_R1:

            for ligne in file_R1:
                if "NO EXONS/GENES PREDICTED IN SEQUENCE" in ligne:
                    return nuc_exon_pred, nuc_intron_pred

                if "Init" in ligne or "Term" in ligne or "Intr" in ligne or "Sngl" in ligne:
                    line = ligne.strip().split(" ")
                    line = list(filter(None, line))
                    start = int(line[3]) - add_nuc
                    stop = int(line[4]) - add_nuc

                    if start >= stop:
                        for nucleotide in range(stop, start + 1):
                            nuc_exon_pred.add(nucleotide)
                    elif start < stop:
                        for nucleotide in range(start, stop + 1):
                            nuc_exon_pred.add(nucleotide)

        end = max(nuc_exon_pred)

        for nuc_intronique in range(-add_nuc + 1, end +add_nuc  + 1):
            if nuc_intronique == 0:
                pass
            else:
                if nuc_intronique not in nuc_exon_pred:
                    nuc_intron_pred.add(nuc_intronique)

        return nuc_exon_pred, nuc_intron_pred
    else:
        print(f"{file} don't exist")
        return nuc_exon_pred, nuc_intron_pred


def create_prediction_geneid(gene, add_nuc, path_pred):
    """
    Generate the structure of the prediction file from Geneid
    :param gene: name of the file .gar
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :param path_pred: path of the prediction files ./Predictions/{soft}/{value}Kb/
    :return:
    """

    soft = "geneid"
    name = gene.replace(".gar", ".fasta")
    gene = f"{name}"
    nuc_intron_pred = set()
    nuc_exon_pred = set()

    if add_nuc != 150:
        add_nuc = int(add_nuc) * 1000
    else:
        add_nuc = 150
    try:
        with open(path_pred + soft + "_" + gene, "r") as file_R1:

            for ligne in file_R1:
                if "exon " in ligne:
                    line = ligne.strip().split("\t")
                    line = list(filter(None, line))
                    tpr = line[0].strip().split(" ")
                    tpr = list(filter(None, tpr))

                    start = int(tpr[1]) - add_nuc
                    stop = int(tpr[2]) - add_nuc

                    if start >= stop:
                        for nucleotide in range(stop, start + 1):
                            nuc_exon_pred.add(nucleotide)
                    elif start < stop:
                        for nucleotide in range(start, stop + 1):
                            nuc_exon_pred.add(nucleotide)

        try:
            end = max(nuc_exon_pred)

            for nuc_intronique in range(-add_nuc, end +add_nuc+ 1):
                if nuc_intronique == 0:
                    pass
                else:
                    if nuc_intronique not in nuc_exon_pred:
                        nuc_intron_pred.add(nuc_intronique)

            return nuc_exon_pred, nuc_intron_pred

        except ValueError:
            return nuc_exon_pred, nuc_intron_pred

    except FileNotFoundError:
        print(f"{gene} not found")
        return nuc_exon_pred, nuc_intron_pred


def create_prediction_glimmer(gene, add_nuc, path_pred):
    """
    Generate the structure of the prediction file from Glimmer
    :param gene: name of the file .gar
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :param path_pred: path of the prediction files ./Predictions/{soft}/{value}Kb/
    :return:
    """

    soft = "glimmer"
    name = gene.replace(".gar", ".fasta")
    gene = f"{name}"

    nuc_exon_pred = set()
    nuc_intron_pred = set()

    if add_nuc != 150:
        add_nuc = int(add_nuc) * 1000
    else:
        add_nuc = 150
    try:
        size = os.path.getsize(path_pred + soft + "_" + name)


        # No predictions
        if size < 100:
            return nuc_exon_pred, nuc_intron_pred
        else:

            with open(path_pred + soft + '_' + gene, "r") as file_R1:

                for ligne in file_R1:
                    if "mRNA\t" in ligne:
                        line = ligne.strip().split("\t")
                        line = list(filter(None, line))

                        start = int(line[3]) - add_nuc
                        stop = int(line[4]) - add_nuc

                        if start >= stop:
                            for nucleotide in range(stop, start + 1):
                                nuc_exon_pred.add(nucleotide)
                        elif start < stop:
                            for nucleotide in range(start, stop + 1):
                                nuc_exon_pred.add(nucleotide)

            end = max(nuc_exon_pred)

            for nuc_intronique in range(-add_nuc, end + add_nuc + 1):
                if nuc_intronique == 0:
                    pass
                else:
                    if nuc_intronique not in nuc_exon_pred:
                        nuc_intron_pred.add(nuc_intronique)

            return nuc_exon_pred, nuc_intron_pred
    except FileNotFoundError:
        print(f"{gene} not found")
        return nuc_exon_pred, nuc_intron_pred


def create_prediction_snap(gene, add_nuc, path_pred):
    """
    Generate the structure of the prediction file from Snap
    :param gene: name of the file .gar
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :param path_pred: path of the prediction files ./Predictions/{soft}/{value}Kb/
    :return:
    """

    soft = "snap"
    name = gene.replace(".gar", ".fasta")
    gene = f"{name}"

    if add_nuc != 150:
        add_nuc = int(add_nuc) * 1000
    else:
        add_nuc = 150

    nuc_exon_pred = set()
    nuc_intron_pred = set()

    try:
        size = os.path.getsize(path_pred + soft + "_" + name)
        # No predictions
        if size < 1:
            return nuc_exon_pred, nuc_intron_pred
        else:
            with open(path_pred + soft + "_" +gene, "r") as file_R1:

                for ligne in file_R1:
                    if "Esngl" in ligne or "Einit" in ligne or "Eterm" in ligne or "Exon" in ligne:
                        line = ligne.strip().split("\t")
                        line = list(filter(None, line))

                        start = int(line[3]) - add_nuc
                        stop = int(line[4]) - add_nuc

                        if start >= stop:
                            for nucleotide in range(stop, start + 1):
                                nuc_exon_pred.add(nucleotide)
                        elif start < stop:
                            for nucleotide in range(start, stop + 1):
                                nuc_exon_pred.add(nucleotide)

            end = max(nuc_exon_pred)

            for nuc_intronique in range(-add_nuc, end +add_nuc + 1):
                if nuc_intronique == 0:
                    pass
                else:
                    if nuc_intronique not in nuc_exon_pred:
                        nuc_intron_pred.add(nuc_intronique)

            return nuc_exon_pred, nuc_intron_pred
    except FileNotFoundError:
        print(f"{gene} not found")
        return nuc_exon_pred, nuc_intron_pred


def comparaison(work_path, gene, soft, add_nuc, path_pred):
    """
    Compare the prediction file with the reference file
    :param work_path: The root of the main project
    :param gene: name of the file .gar
    :param soft: Name of the predictor
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :param path_pred: path of the prediction files ./Predictions/{soft}/{value}Kb/
    :return:
    """
    if soft == "augustus":
        exon_pred_exon = create_prediction_augustus(gene, add_nuc, path_pred)[0]
        intron_pred_exon = create_prediction_augustus(gene, add_nuc, path_pred)[2]
        exon_pred_cds = create_prediction_augustus(gene, add_nuc, path_pred)[1]
        intron_pred_cds = create_prediction_augustus(gene, add_nuc, path_pred)[3]

        if add_nuc == 150:
            output = work_path + f"Resultats/nucleotide/{soft}/150bp/"
        else:
            output = work_path + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/"

        os.makedirs(output, exist_ok=True)

    elif soft == "genscan":
        exon_predit = create_prediction_genscan(gene, add_nuc, path_pred)[0]
        intron_predit = create_prediction_genscan(gene, add_nuc, path_pred)[1]

        if add_nuc == 150:
            output = work_path + f"Resultats/nucleotide/{soft}/150bp/"
        else:
            output = work_path + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/"
        os.makedirs(output, exist_ok=True)

    elif soft == "geneid":
        exon_predit = create_prediction_geneid(gene, add_nuc, path_pred)[0]
        intron_predit = create_prediction_geneid(gene, add_nuc, path_pred)[1]

        if add_nuc == 150:
            output = work_path + f"Resultats/nucleotide/{soft}/150bp/"
        else:
            output = work_path + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/"
        os.makedirs(output, exist_ok=True)

    elif soft == "glimmer":
        exon_predit = create_prediction_glimmer(gene, add_nuc, path_pred)[0]
        intron_predit = create_prediction_glimmer(gene, add_nuc, path_pred)[1]

        if add_nuc == 150:
            output = work_path + f"Resultats/nucleotide/{soft}/150bp/"
        else:
            output = work_path + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/"
        os.makedirs(output, exist_ok=True)

    elif soft == "snap":
        exon_predit = create_prediction_snap(gene, add_nuc, path_pred)[0]
        intron_predit = create_prediction_snap(gene, add_nuc, path_pred)[1]

        if add_nuc == 150:
            output = work_path + f"Resultats/nucleotide/{soft}/150bp/"
        else:
            output = work_path + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/"
        os.makedirs(output, exist_ok=True)

    else:
        print("ERROR in program name, please choose : augustus, genscan, geneid, glimmer or snap")
        exit()

    filename = gene.replace(".gar", "")

    # Calculate the different metrics
    for transcrit, nuc_ref in create_reference(work_path, gene, soft).items():
        if soft == "augustus":
            tp_exon = set()
            fp_exon = set()
            tn_exon = set()
            fn_exon = set()

            tp_cds = set()
            fp_cds = set()
            tn_cds = set()
            fn_cds = set()

            # EXON
            for nuc_exon in nuc_ref[0]:
                # TP
                if nuc_exon in exon_pred_exon:
                    tp_exon.add(nuc_exon)
                # FN
                elif nuc_exon not in exon_pred_exon:
                    fn_exon.add(nuc_exon)

            for nuc_intron in nuc_ref[1]:
                # TN
                if nuc_intron in intron_pred_exon:
                    tn_exon.add(nuc_intron)
                # FP
                elif nuc_intron not in intron_pred_exon:
                    fp_exon.add(nuc_intron)

            # CDS
            for nuc_exon in nuc_ref[0]:
                # TP
                if nuc_exon in exon_pred_cds:
                    tp_cds.add(nuc_exon)
                # FN
                elif nuc_exon not in exon_pred_cds:
                    fn_cds.add(nuc_exon)

            for nuc_intron in nuc_ref[1]:
                # TN
                if nuc_intron in intron_pred_cds:
                    tn_cds.add(nuc_intron)
                # FP
                elif nuc_intron not in intron_pred_cds:
                    fp_cds.add(nuc_intron)

            tp_exon = len(tp_exon)
            fp_exon = len(fp_exon)
            tn_exon = len(tn_exon)
            fn_exon = len(fn_exon)
            tp_cds = len(tp_cds)
            fp_cds = len(fp_cds)
            tn_cds = len(tn_cds)
            fn_cds = len(fn_cds)


            # calculate the sensitivity
            sn_exon = sensibility(tp_exon, fn_exon)
            # calculate the specificity
            sp_exon = specificity(tp_exon, fp_exon)
            # calculate the F1 score
            f1_exon = f1_score(tp_exon, fp_exon, fn_exon)

            # calculate the sensitivity
            sn_cds = sensibility(tp_cds, fn_cds)
            # calculate the specificity
            sp_cds = specificity(tp_cds, fp_cds)
            # calculate the F1 score
            f1_cds = f1_score(tp_cds, fp_cds, fn_cds)

            if f1_exon > f1_cds:
                # Writting file
                with open(output + soft + "_all_metrics.csv", "a") as file_W1:
                    file_W1.write(
                        f"{filename.replace('.fasta', '')}_exon;{str(sn_exon)};{str(sp_exon)};{str(f1_exon)};{str(tp_exon)};{str(fp_exon)};{str(tn_exon)};{str(fn_exon)}\n")
            else:
                # Writting file
                with open(output + soft + "_all_metrics.csv", "a") as file_W1:
                    file_W1.write(
                        f"{filename.replace('.fasta', '')}_cds;{str(sn_cds)};{str(sp_cds)};{str(f1_cds)};{str(tp_cds)};{str(fp_cds)};{str(tn_cds)};{str(fn_cds)}\n")

        # Other soft
        else:
            tp = set()
            fp = set()
            tn = set()
            fn = set()

            for nuc_exon in nuc_ref[0]:
                # TP
                if nuc_exon in exon_predit:
                    tp.add(nuc_exon)
                # FN
                elif nuc_exon not in exon_predit:
                    fn.add(nuc_exon)

            for nuc_intron in nuc_ref[1]:
                # TN
                if nuc_intron in intron_predit:
                    tn.add(nuc_intron)
                # FP
                elif nuc_intron not in intron_predit:
                    fp.add(nuc_intron)

            tp = len(tp)
            fp = len(fp)
            tn = len(tn)
            fn = len(fn)

            # calculate the sensitivity
            sn = sensibility(tp, fn)
            # calculate the specificity
            sp = specificity(tp, fp)
            # calculate the F1 score
            f1 = f1_score(tp, fp, fn)

            # Writing files
            with open(output + soft + "_all_metrics.csv", "a") as file_W1:
                file_W1.write(
                    f"{filename.replace('.fasta', '')};{str(sn)};{str(sp)};{str(f1)};{str(tp)};{str(fp)};{str(tn)};{str(fn)}\n")


def create_header(work_path, output, soft, add_nuc):
    """
    Create an header in the different output file
    :param work_path: The root of the main project
    :param output: path of the output
    :param soft: Name of the predictor
    :param add_nuc: Number of nucleotides flanked upstream and downstream of gene
    :return:
    """
    if add_nuc == 150:
        output_dir = work_path + f"Resultats/nucleotide/{soft}/150bp/"
    else:
        output_dir = work_path + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/"

    os.makedirs(output_dir, exist_ok=True)

    with open(output, "w") as file_W1:
        file_W1.write("Filename;Sensibility;Specificity;F1_score;TP;FP;TN;FN\n")


def create_header_final(work_path, output):
    """
    Create an header in the final output file
    :param work_path: The root of the main project
    :param output: path of the output
    :return:
    """
    output_dir = work_path + f"Resultats/nucleotide/"
    os.makedirs(output_dir, exist_ok=True)
    with open(output, "a") as file_W2:
        file_W2.write(f"SN;SP;F1;\n")


def calcule_moyenne(input_f, output):
    """
    Calculates the average for all the metrics
    :param input_f: {program}_all_metrics.csv
    :param output: {program}_results.csv
    :return:
    """

    mean_sensi = []
    mean_speci = []
    mean_f1 = []

    with open(input_f, "r") as file_R1:
        for i, ligne in enumerate(file_R1):

            if i == 0:
                pass
            else:
                ligne = ligne.strip().split(";")

                sensi = float(ligne[1])
                speci = float(ligne[2])
                f1 = float(ligne[3])

                mean_sensi.append(sensi)
                mean_speci.append(speci)
                mean_f1.append(f1)

    msensi = str(round(moyenne(mean_sensi), 2))
    mspeci = str(round(moyenne(mean_speci), 2))
    mf1 = str(round(moyenne(mean_f1), 2))

    with open(output, "a") as file_W1:
        file_W1.write(f"{msensi};{mspeci};{mf1}\n")



def main(workpath, soft, add_nuc):

    #  Creation of paths for output files.
    if add_nuc == 150:
        output_all_metrics = workpath + f"Resultats/nucleotide/{soft}/150bp/{soft}_all_metrics.csv"
        output_final = workpath + f"Resultats/nucleotide/{soft}/150bp/{soft}_results.csv"
        path_pred = workpath + f"Predictions/{soft}/150bp/"

    else:
        output_all_metrics = workpath + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/{soft}_all_metrics.csv"
        output_final = workpath + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/{soft}_results.csv"
        path_pred = workpath + f"Predictions/{soft}/{str(add_nuc)}Kb/"

    # Creation of header
    create_header(workpath, output_all_metrics, soft, add_nuc)
    create_header_final(workpath, output_final)

    # Generation of gene list
    liste_file = liste_gene(workpath, soft, add_nuc)

    # Launch of comparisons
    for gene in tqdm(liste_file):
        # print(gene)
        comparaison(workpath, gene, soft, add_nuc, path_pred)

    if add_nuc == 150:
        calcule_moyenne(workpath + f"Resultats/nucleotide/{soft}/150bp/{soft}_all_metrics.csv", workpath + f"Resultats/nucleotide/{soft}/150bp/{soft}_results.csv")
    else:
        calcule_moyenne(workpath  + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/{soft}_all_metrics.csv", workpath + f"Resultats/nucleotide/{soft}/{str(add_nuc)}Kb/{soft}_results.csv")


if __name__ == '__main__':
    # workpath is the root directory where you download the git repository
    workpath = "my_path/Benchmark_study/"
    l_soft=["augustus", "genscan","geneid", "glimmer","snap"]
    l_nuc = [150, 2,4,6,8,10]

    for soft in l_soft:
        for add_nuc in l_nuc:
            main(workpath, soft, add_nuc)





