import torch
import cobra
import datetime
import pandas as pd
import numpy as np
import subprocess
import re
import sys
import os
import math
import shutil
sys.path.append(r'./script/')
from ECMpy_function import *
#from ECMpy_function import get_enzyme_constraint_model



def Model_preview(model_file):
    bigg_met_file = './data/bigg_models_metabolites.txt'
    result = Determine_suitable_ecGEM(model_file, bigg_met_file)
    return result


def get_kcat_via_DLKcat(model_file, kcat_analysis_folder):

    create_file(dlkcat_folder)
    # input files
    sbml_path = model_file

    # output files
    gene_subnum_path = "%sgene_subnum.csv" % dlkcat_folder
    sub_description_path = '%sget_gene_subunitDescription.csv' % dlkcat_folder
    inchikey_list_file = '%sinchikey_list.csv' % dlkcat_folder
    inchikey_list_smilesfile = '%sinchikey_list_smiles.csv' % dlkcat_folder
    comdf_file = '%scomdf.csv' % dlkcat_folder
    DLouputdf_file = '%sDLoutput.tsv' % dlkcat_folder
    metdf_outfile = '%smetabolites_reactions_gpr_similes_prosequence_mass_dropna.csv' % dlkcat_folder
    metabolites_reactions_gpr_file = '%smetabolites_reactions_gpr.csv' % dlkcat_folder
    prodf_file = '%sprodf.csv' % dlkcat_folder
    DLinput_file = '%sDLinput.tsv' % dlkcat_folder
    DL_reaction_kact_mw_file = '%sreaction_kcat_MW.csv' % dlkcat_folder

    model = cobra.io.read_sbml_model(sbml_path)

    #NB: the following steps can be run in one go, or one by one.
    #Step 1: subunit number of each reaction,
    # Run this part only once --> then comment out:
    # starttime = datetime.datetime.now()
    # print("Starting to fetch subunit number of each enzyme")
    # get_gene_subunitDescription(sub_description_path, model)  # Download from the UniProt API, run it once.
    # subbnumdf = get_subunit_number(sub_description_path, gene_subnum_path)
    # print("Calculation done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Step 2: Convert metabolites bigg id to smiles
    # starttime = datetime.datetime.now()
    # # Step 2: convert metbolites bigg id to smiles
    # print("Starting to convert metbolites bigg id to smiles...")
    # metdf_name = get_met_bigg_id(model)
    # inchkeydf = convert_bigg_met_to_inchikey(metdf_name['met'], inchikey_list_file)  # from BIGG
    # # inchkeydf = pd.read_csv('./data/inchikey_list.csv')
    # smilesdf = convert_inchikey_to_smiles(inchkeydf, inchikey_list_smilesfile)  # from pubchem
    # print("Converting done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Step 3: get protein sequence and mass in model
    # starttime = datetime.datetime.now()
    # # Step 3: get protein sequence and mass in model
    # print("Starting to get protein sequence and mass in model...")
    # subbnumdf = pd.read_csv(gene_subnum_path)
    # prodf = get_model_protein_sequence_and_mass(model, subbnumdf, prodf_file)
    # print("Getting done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Quick and dirty solution --> update the prodf_file so that subunitmass matches the mass column
    # prodf = pd.read_csv(prodf_file)
    # prodf["subunitmass"] = prodf["subunitmass"].fillna(prodf["mass"])
    # print(prodf.head())
    # upd_prodf_file = '%sprodf_upd.csv' % dlkcat_folder
    # NB: I replaced the old prodf file, with the updated one, so this code cannot be rerun.
    # prodf.to_csv(prodf_file, index=False)

    #Step 4: split the substrate of reactions to match the gene
    # starttime = datetime.datetime.now()
    # print("Starting to split the substrate of reactions to match the gene...")
    # spdf = split_substrate_to_match_gene(model, metabolites_reactions_gpr_file)
    # print("Splitting done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)
    # exit()

    ## Step 5: combine the reaction--substrate--gene--protein_sequnce--mass and formate DLKcat input file
    # starttime = datetime.datetime.now()
    # print("Starting to combine data...")
    # metdf_name = get_met_bigg_id(model)
    # smilesdf = pd.read_csv(inchikey_list_smilesfile)
    # spdf = pd.read_csv(metabolites_reactions_gpr_file)
    # prodf = pd.read_csv(prodf_file)
    # comdf = combine_reactions_smiles_sequence(spdf, smilesdf, prodf, comdf_file)
    # DLinputdf = generate_DLKCAT_input(comdf, metdf_name, metdf_outfile, DLinput_file)
    # print("Combining done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    #Step 6: use DLKcat calculate kcat
    # starttime = datetime.datetime.now()
    # print("Starting to Use DLKcat calculate kcat...")
    # cmd_str = f"python ./script/prediction_for_input.py ./ {dlkcat_folder}DLinput.tsv {dlkcat_folder}DLoutput.tsv"
    # subprocess.run(cmd_str, shell=True)
    # print("DLKcat done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)
    # Step 7: get the kcat_mw file
    starttime = datetime.datetime.now()

    print("Starting to get reaction kcat_mw for model......")
    DLouputdf = pd.read_csv(DLouputdf_file, sep='\t')
    comdf = pd.read_csv(comdf_file)
    DL_reaction_kact_mw = DL_kcat_mw_calculation(DLouputdf, comdf)
    DL_reaction_kact_mw.to_csv(DL_reaction_kact_mw_file, index=False)

    endtime = datetime.datetime.now()
    print(endtime - starttime)

    # everything has finished without errors (biological interpretation needs to be checked)
    return None


def get_kcat_via_autopacman(model_file, autopacmen_folder):
    kcat_gap_fill = 'mean'  # 'mean'#'median'
    reaction_gap_fill = 'mean'
    sbml_path = model_file
    organism = "Streptomyces clavuligerus"
    project_name = "iDG1237_%s" % kcat_gap_fill
    create_file(autopacmen_folder)
    protein_kcat_database_path = "none"
    bigg_metabolites_file = "./data/bigg_models_metabolites.txt"  # date:20230629 http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt
    brenda_textfile_path = "../input_files/brenda_2026_1.txt"  # date:20260318 https://www.brenda-enzymes.org/download.php
    uniprot_data_file = './data/uniprot_data_accession_key.json'  # from uniprot

    # output files
    brenda_json_path = "%skcat_database_brenda.json" % autopacmen_folder
    brenda_json_path2 = "%ssa_database_brenda.json" % autopacmen_folder
    sabio_rk_json_path = "%skcat_database_sabio_rk.json" % autopacmen_folder
    bigg_id_name_mapping_path = "%sbigg_id_name_mapping.json" % autopacmen_folder
    brenda_output_json_path = "%skcat_database_brenda_for_model.json" % autopacmen_folder
    combined_output_path = "%skcat_database_combined.json" % autopacmen_folder
    sub_description_path = '%sget_gene_subunitDescription.csv' % autopacmen_folder
    gene_subnum_path = "%sgene_subnum.csv" % autopacmen_folder
    reaction_mw_path = "%sreaction_mw.json" % autopacmen_folder
    reaction_kcat_mw_path = '%sreaction_kcat_MW.csv' % autopacmen_folder

    # step 1: get bigg metabolite
    # starttime = datetime.datetime.now()
    # print("Starting to deal BIGG metabolites text file...")
    # parse_bigg_metabolites_file(bigg_metabolites_file, autopacmen_folder)
    # print("BIGG metabolites text file done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # step 2: BRENDA kcat
    # starttime = datetime.datetime.now()
    # print("Starting to deal BRENDA textfile...")
    # parse_brenda_textfile(brenda_textfile_path, autopacmen_folder, brenda_json_path, brenda_json_path2)
    # print("BRENDA textfile done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Step 3: Select Brenda kcat for model
    # starttime = datetime.datetime.now()
    # print("Starting to deal brenda json for model...")
    # parse_brenda_json_for_model(sbml_path, brenda_json_path, brenda_output_json_path)
    # print("BRENDA json for model done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Step 4: SABIO-RK kcat for model
    # starttime = datetime.datetime.now()
    # print("Starting EC numbers kcat search in SABIO-RK...")
    # parse_sabio_rk_for_model_with_sbml(sbml_path, sabio_rk_json_path, bigg_id_name_mapping_path)
    # print("SABIO-RK done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Step 5: Brenda and SABIO-RK kcat combined
    # starttime = datetime.datetime.now()
    # print("Combining kcat database...")
    # create_combined_kcat_database(sabio_rk_json_path, brenda_output_json_path, combined_output_path)
    # print("Combining kcat database done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Step 6: subunit number of each reaction
    # starttime = datetime.datetime.now()
    # print("Starting to fetch subunit number of each enzyme")
    # if re.search('\.xml', sbml_path):
    #     model = cobra.io.read_sbml_model(sbml_path)
    # elif re.search('\.json', sbml_path):
    #     model = cobra.io.json.load_json_model(sbml_path)
    # get_gene_subunitDescription(sub_description_path, model)  # 从uniprot的api下载，运行一次就行
    # subbnumdf = get_subunit_number(sub_description_path, gene_subnum_path)
    # print("Calculation done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Step 7: get mw for model gene (must be uniprot ID)
    # starttime = datetime.datetime.now()
    # print("Starting UniProt ID<->Protein mass search using UniProt...")
    # # get_protein_mass_mapping_from_local(sbml_path, autopacmen_folder, project_name, uniprot_data_file)
    # get_protein_mass_mapping_with_sbml(sbml_path, autopacmen_folder, project_name)
    #
    # get_reaction_mw(sbml_path, autopacmen_folder, project_name, reaction_mw_path, gene_subnum_path)
    # print("Protein ID<->Mass mapping done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Step 8: kcat assignment for model(include sa)
    # starttime = datetime.datetime.now()
    # print("Starting to assign kcat for model...")
    # get_reactions_kcat_mapping(sbml_path, autopacmen_folder, project_name, organism, combined_output_path,
    #                            brenda_json_path2, reaction_mw_path, protein_kcat_database_path, kcat_gap_fill)
    # print("kcat assignment done!")
    # print()
    #
    # endtime = datetime.datetime.now()
    # print(endtime - starttime)

    # Step 9: get_reaction_kcat_mw for model
    starttime = datetime.datetime.now()
    print("Starting to get reaction kcat_mw for model...")
    if re.search('\.xml', sbml_path):
        model = cobra.io.read_sbml_model(sbml_path)
    elif re.search('\.json', sbml_path):
        model = cobra.io.json.load_json_model(sbml_path)
    get_reaction_kcat_mw(model, autopacmen_folder, project_name, reaction_gap_fill, gene_subnum_path,
                         reaction_kcat_mw_path)
    print("Reaction kcat_mw done!")

    endtime = datetime.datetime.now()
    print(endtime - starttime)

    return None


def get_ecModel_using_ECMpy(input_file, analysis_folder, additional_kcat_values=None):
    sbml_path = input_file
    taxonom_id = 1901
    kcat_folder = analysis_folder
    reaction_kcat_MW_file = f"{kcat_folder}reaction_kcat_MW.csv"

    # this part is recommemend but Sclav is not part of the pax-db database.
    # gene_abundance_colname = 'abundance'
    # # gene_abundance_file = './data/gene_abundance.csv'  # downolad from https://pax-db.org/download
    # # # The enzyme mass fraction,such as 0.405
    # # # f=calculate_f_v2(sbml_path, gene_abundance_file,gene_abundance_colname,taxonom_id)

    f = 0.405
    # Initial parameters
    ptot = 0.56  # The total protein fraction in cell.
    sigma = 1  # The approximated saturation of enzyme.e.g.,0.5/1.
    lowerbound = 0  # Lowerbound  of enzyme concentration constraint.
    upperbound = round(ptot * f * sigma, 3)  # total enzyme
    ecModel_output_file = "./model/iDG1237_updated_irr_enz_constraint_v2.json"

    #NB: the result of f*ptot is estimated now, but can later still be altered, also in the model, so this function
    # does not have to be rerun to update the size of the proteome pool

    # Get ecModel
    #only do this part once, then comment out
    # trans_model2enz_json_model_split_isoenzyme(sbml_path, reaction_kcat_MW_file, f, ptot, sigma, lowerbound, upperbound,
    #                                            ecModel_output_file)

    # the resulting model still contains NaN values. They need be dealt with. They should be simply put to zero, because
    # then their usage would be free. So: use the file from Sepideh with Brenda and GECKO values. Put remaining values
    # at 65.9.
    ecModel_output_file_upd = "./model/iDG1237_updated_irr_enz_constraint_v4.json"

    kcat_df = pd.read_csv(additional_kcat_values, sep='\t')
    mw_df = pd.read_csv(reaction_kcat_MW_file)
    enz_model = get_enzyme_constraint_model(ecModel_output_file)
    counter = 0
    for reac in enz_model.reactions:
        upd_id = reac.id
        if "_reverse" in upd_id:
            upd_id = upd_id.replace("reverse", "REV")
        if '_num' in upd_id:
            upd_id = upd_id.replace('_num', '_EXP_')
        if isinstance(reac.kcat, (float, np.floating)) and np.isnan(reac.kcat):
            if upd_id in kcat_df["rxns"].values:
                try:
                    kcat = kcat_df.loc[kcat_df["rxns"] == upd_id, "kcats"].values[0]
                except IndexError:
                    print(upd_id)
                counter += 1
            else:
                if 'transport' in reac.name:
                    kcat = ""
                else:
                    kcat = 65.9
            # calculate kcat_MW: DLoutputdf_rex['kcat_mw'] = DLoutputdf_rex['Kcat value (1/s)'] * 3600 * 1000 / DLoutputdf_rex['totalmass']
            reac.kcat = kcat
            if kcat == "":
                reac.kcat_MW = ""
            else:
                MW = mw_df.loc[mw_df["reactions"] == reac.id, "MW"].values[0]
                if isinstance(MW, (float, np.floating)) and np.isnan(MW):
                    reac.kcat_MW = ""
                    #print(f"MW of {upd_id} is Nan")
                else:
                    kcat_mw = kcat * 3600 * 1000 / MW
                    reac.kcat_MW = kcat_mw
        else:
            if upd_id in kcat_df["rxns"].values:
                try:
                    kcat = kcat_df.loc[(kcat_df["rxns"] == upd_id) & (kcat_df["kcatSource"] == "brenda"), "kcats"].values[0]
                except IndexError:
                    continue
                try:
                    MW = mw_df.loc[mw_df["reactions"] == reac.id, "MW"].values[0]
                except IndexError:
                    print(f"MW problem for {reac.id}")
                    continue
                if isinstance(MW, (float, np.floating)) and np.isnan(MW):
                    reac.kcat_MW = ""
                    # print(f"MW of {upd_id} is Nan")
                else:
                    reac.kcat = kcat
                    kcat_mw = kcat * 3600 * 1000 / MW
                    reac.kcat_MW = kcat_mw
                print(f"{kcat_mw} is the new value for {reac.id}")

    # use the old model as the frame-work and fill it in:
    dictionary_model = json_load(ecModel_output_file)
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        dictionary_model['reactions'][eachreaction]['kcat'] = enz_model.reactions.get_by_id(reaction_id).kcat
        dictionary_model['reactions'][eachreaction]['kcat_MW'] = enz_model.reactions.get_by_id(reaction_id).kcat_MW
        kcat_mw = enz_model.reactions.get_by_id(reaction_id).kcat_MW
        if isinstance(kcat_mw, (float, np.floating)) and np.isnan(kcat_mw):
            print(reaction_id)
    json_write(ecModel_output_file_upd, dictionary_model)
    return None


def get_ecModel_for_autopacman(input_file, analysis_folder):
    sbml_path = input_file
    taxonom_id = 1901
    kcat_folder = analysis_folder
    reaction_kcat_MW_file = f"{kcat_folder}reaction_kcat_MW.csv"

    # this part is recommemend but Sclav is not part of the pax-db database.
    # gene_abundance_colname = 'abundance'
    # # gene_abundance_file = './data/gene_abundance.csv'  # downolad from https://pax-db.org/download
    # # # The enzyme mass fraction,such as 0.405
    # # # f=calculate_f_v2(sbml_path, gene_abundance_file,gene_abundance_colname,taxonom_id)

    f = 0.405
    # Initial parameters
    ptot = 0.56  # The total protein fraction in cell.
    sigma = 1  # The approximated saturation of enzyme.e.g.,0.5/1.
    lowerbound = 0  # Lowerbound  of enzyme concentration constraint.
    upperbound = round(ptot * f * sigma, 3)  # total enzyme
    ecModel_output_file = "./model/iDG1237_updated_irr_enz_constraint_autopacman_v1.json"

    # Get ecModel
    #only do this part once, then comment out
    trans_model2enz_json_model_split_isoenzyme(sbml_path, reaction_kcat_MW_file, f, ptot, sigma, lowerbound, upperbound,
                                               ecModel_output_file)
    return None


if __name__ == '__main__':
    #input_model_file = "../Model_files/NewModel_CleanedFormula_Hbalance3.xml"
    #input_model_file = "../Model_files/iDG1237.xml"
    input_model_file = "../Model_files/iDG1237_updated_with_uniprot.xml"

    #step zero: do the model check for Emcpy
    # model_check = Model_preview(input_model_file)
    # print(model_check)
    #Output is: yes, it is suitable

    #Step 1: DLKcat
    dlkcat_folder = "./analysis/get_kcat_mw_by_DLKcat_v1_iDG1237/"
    #get_kcat_via_DLKcat(input_model_file, dlkcat_folder)


    database_kcat_values = "../input_files/kcatList_merged_GECKO.tsv"
    # Step 2: get ecModel with ECMpy
    #get_ecModel_using_ECMpy(input_model_file, dlkcat_folder, additional_kcat_values=database_kcat_values)

    # Alternative to step 1: use Autopacman, to see what kcat values we get from that.
    out_folder = "./analysis/get_kcat_mw_by_AutoPACMEN_v1_iDG1237/"
    #get_kcat_via_autopacman(input_model_file, out_folder)

    #updated step 2: get ecModel with ECMpy
    get_ecModel_for_autopacman(input_model_file, out_folder)



