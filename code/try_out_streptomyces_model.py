import os.path

import cobra
import pandas as pd

def compare_models(ori_model_file, new_model_file):
    ori_model = cobra.io.read_sbml_model(ori_model_file)
    new_model = cobra.io.read_sbml_model(new_model_file)
    # print(ori_model.objective.expression)
    # print(new_model.objective.expression)
    ori_sol = ori_model.optimize()
    new_sol = new_model.optimize()
    ori_medium = ori_model.medium
    # for nutrient in ori_medium:
    #     print(nutrient, ori_sol.fluxes[nutrient])
    new_medium = new_model.medium
    # for nutrient in new_medium:
    #     print(nutrient, new_sol.fluxes[nutrient])
    #print(dir(new_model.reactions.EX_SPC_e))
    # print(new_model.reactions.EX_SPC_e.name)
    # print(new_model.metabolites.SP_c.reactions)
    # #print(new_model.metabolites.SP_c.name)
    #
    # print(new_model.reactions.SOYPEPTtr.name)
    # print(new_model.reactions.SOYPEPTtr.reaction)
    # print(new_model.reactions.SOYPEPT.name)
    # print(new_model.reactions.SOYPEPT.reaction)

    # print(ori_model.reactions.ATPM.bounds)
    # print(new_model.reactions.ATPM.bounds)

    # print(ori_model.reactions.Growth.reaction)
    # print(new_model.reactions.Growth.reaction)

    new_model.reactions.EX_pi_e.lower_bound = -0.14
    #new_model.reactions.EX_nh4_e.upper_bound = 1000
    new_model.reactions.EX_nh4_e.upper_bound = 0
    # new_sol = new_model.optimize()
    # for nutrient in new_medium:
    #     print(nutrient, new_sol.fluxes[nutrient])
    # print('EX_pi_e', new_sol.fluxes["EX_pi_e"])
    # print(new_sol)
    #
    # print(new_model.reactions.EX_o2_e.bounds)

    # for metabolite in ori_model.metabolites:
    #     meta_name = ori_model.metabolites.get_by_id(metabolite.id).name
    #     if "clavu" in meta_name.lower():
    #         print(metabolite.id, meta_name)

    # print(ori_model.metabolites.get_by_id("clava[e]").reactions)
    # print(ori_sol.fluxes['EX_clav_e'])
    # print(new_model.metabolites.get_by_id("clava_e").reactions)
    # print(new_sol.fluxes['EX_clav_e'])

    #change the objective to ATPM
    new_model.objective = "ATPM"
    new_model.reactions.EX_pi_e.lower_bound = 0

    upd_sol = new_model.optimize()
    print(upd_sol)
    print(upd_sol.fluxes['EX_glyc_e'])
    print(upd_sol.fluxes['Growth'])
    #new_model.reactions.EX_pi_e.lower_bound = -0.14

    #check what happens when all medium bounds are closed
    # for reac in new_model.reactions:
    #     if reac.id.startswith('EX_'):
    #         reac.lower_bound = 0
    # upd_sol = new_model.optimize()
    # print(upd_sol)

    # for meta in ori_model.metabolites:
    #     if ";" in meta.formula:
    #         print(meta.id)
    #         print(meta.formula)

    # print(ori_model.metabolites.get_by_id("ocdcea[c]"))
    # print(ori_model.metabolites.get_by_id("ocdcea[e]"))
    # print(ori_model.metabolites.get_by_id("ocdca[c]"))
    # print(ori_model.metabolites.get_by_id("ocdca[e]"))

    # print(new_model.metabolites.get_by_id("ocdcea_c").name)
    # print(new_model.metabolites.get_by_id("ocdcea_c").formula)
    # print(new_model.metabolites.get_by_id("ocdcea_e"))
    # print(new_model.metabolites.get_by_id("ocdcea_p"))

    # print(dir(new_model.reactions.get_by_id("CITL")))
    # print(new_model.reactions.get_by_id("CITL").metabolites)
    #
    # for ori_reaction in ori_model.reactions:
    #     try:
    #         new_reaction = new_model.reactions.get_by_id(ori_reaction.id)
    #     except KeyError:
    #         print(ori_reaction.id)
    #         continue
    # ori_gene_name_list = []
    # for gene in ori_model.genes:
    #     gene_name = gene.notes
    #     ori_gene_name_list.append(gene_name)
    # #
    # new_gene_name_list = []
    # for gene in new_model.genes:
    #     gene_name = gene.notes
    #     new_gene_name_list.append(gene_name)
    #
    # print(len(ori_gene_name_list))
    # print(len(new_gene_name_list))
    #
    # print(ori_gene_name_list[:5])
    # print(new_gene_name_list[:5])
    ori_meta_list = []
    for meta in ori_model.metabolites:
        ori_meta_list.append(meta.id)
    new_meta_list = []
    for meta in new_model.metabolites:
        if '_c' in meta.id or '_p' in meta.id or '_e' in meta.id:
            try:
                metabolite, compartment = meta.id.rsplit('_', 1)
                new_id = f"{metabolite}[{compartment}]"
                new_meta_list.append(new_id)
            except ValueError:
                print(meta.id)
        else:
            new_meta_list.append(meta.id)

    print(len(ori_meta_list))
    print(len(new_meta_list))

    first_list = [meta for meta in ori_meta_list if meta not in new_meta_list]
    second_list = [meta for meta in new_meta_list if meta not in ori_meta_list]

    print("meta that are in ori model, but not in new model")
    for elem in first_list:
        print(elem)

    print("meta that are in new model, but not in ori model")

    for elem in second_list:
        print(elem)

    df_rows = []
    # save all metabolite information
    for meta_id in second_list:
        id_elem = meta_id.split('[')
        upd_id = id_elem[0] + '_' + id_elem[1][0]
        meta = new_model.metabolites.get_by_id(upd_id)
        df_rows.append({
            "id": meta.id,
            "name": meta.name,
            "formula": meta.formula,
            "charge": meta.charge,
            "compartment": meta.compartment,
            "annotation": meta.annotation,
            "notes": meta.notes,
        })
    df = pd.DataFrame(df_rows)
    print(df)
    if not os.path.exists("SclavGEM/data/data_to_update_published_model/new_metabolites.tsv"):
        df.to_csv("SclavGEM/data/data_to_update_published_model/new_metabolites.tsv", sep="\t", index=False)

    # count = 0
    # for meta in ori_model.metabolites:
    #     print(meta.formula)
    #     if meta.formula == "":
    #         count += 1
    # print(count)

    return None


def try_out_clav_production_cleaned_formula_model(new_model_file):
    new_model = cobra.io.read_sbml_model(new_model_file)
    # current result: we have CLAV production, both for obj growth and obj ATPM
    new_model.objective = "ATPM"
    #new_model.reactions.EX_pi_e.lower_bound = -0.14
    #new_model.reactions.EX_pi_e.lower_bound = -2
    new_model.reactions.EX_o2_e.bounds = -4, -4
    new_model.reactions.EX_SPC_e.lower_bound = 0
    new_model.reactions.EX_nh4_e.lower_bound = -1000
    sol = new_model.optimize()
    print(sol)
    print(sol.fluxes['EX_o2_e'])
    print(sol.fluxes['EX_glyc_e'])
    #print('ADK flux', print(sol.fluxes['ADK2']))
    #print(new_model.metabolites.atp_c.summary())
    print(new_model.metabolites.o2_c.summary())

    # limit oxygen consuming complexes
    new_model.reactions.CYTBD.upper_bound = 1
    new_model.reactions.CYTBD2pp.upper_bound = 1
    new_model.reactions.CYTBD2.upper_bound = 1
    new_model.reactions.CYTBDpp_1.upper_bound = 1
    new_model.reactions.CYTB_B2.upper_bound = 1
    new_model.reactions.CYO1b.upper_bound = 1
    new_model.reactions.DESAT16.upper_bound = 0
    new_model.reactions.DESAT18.upper_bound = 0
    new_model.reactions.ROO.upper_bound = 0
    new_model.reactions.NOX.upper_bound = 0
    new_model.reactions.GLYCTO1.upper_bound = 0
    new_model.reactions.ASPO1.upper_bound = 0
    new_model.reactions.DALAOX.upper_bound = 0
    new_model.reactions.PYAM5PO.upper_bound = 0
    new_model.reactions.MEHLER.upper_bound = 0
    #after fixing o2 bounds:
    new_model.reactions.FEROpp.upper_bound = 0
    new_model.reactions.GGPTRCO.upper_bound = 0
    #after moving to minimal medium without soy protein
    new_model.reactions.HGNTOR.upper_bound = 0
    new_model.reactions.get_by_id("34HPPOR").upper_bound = 0
    new_model.reactions.get_by_id("FAS181").upper_bound = 0
    new_model.reactions.get_by_id("FAS161").upper_bound = 0
    new_model.reactions.get_by_id("FDMO2").upper_bound = 0
    new_model.reactions.get_by_id("FDMO_1").upper_bound = 0
    new_model.reactions.get_by_id("LYSMO").upper_bound = 0
    new_model.reactions.get_by_id("INOSTO").upper_bound = 0
    new_model.reactions.get_by_id("PYROX").upper_bound = 0
    sol = new_model.optimize()
    print('results after modulating O2-consuming reactions')
    print(sol)
    print(sol.fluxes['EX_o2_e'])
    print(sol.fluxes['EX_glyc_e'])
    print("clav flux is: " , sol.fluxes['EX_clav_e'])
    print(sol.fluxes['EX_pi_e'])
    # print(new_model.metabolites.o2_c.summary())
    # print(new_model.metabolites.o2_p.summary())
    # print(new_model.metabolites.clavam_c.summary())
    # print(new_model.metabolites.clavad_c.summary())
    # print(new_model.metabolites.cmeclav_c.summary())
    # print(new_model.metabolites.fomclav_c.summary())
    # print(new_model.metabolites.hmclav_c.summary())
    # print(new_model.metabolites.fclav_c.summary())
    # print(new_model.metabolites.clavcar_c.summary())
    print("clavam-2-carboxylate flux is: ", sol.fluxes['EX_clavcar_e'])

    #print(new_model.metabolites.atp_c.summary(solution=sol))

    for reac in new_model.reactions:
        if reac.id.startswith('EX_'):
            flux = sol.fluxes[reac.id]
            if flux != 0:
                print(reac.id, flux)
    #print(new_model.metabolites.h_c.summary(solution=sol))
    print(new_model.metabolites.o2_c.summary(solution=sol))
    return None


def try_out_clav_production_updated_uniprot_model(new_model_file):
    new_model = cobra.io.read_sbml_model(new_model_file)
    # current result: we have CLAV production, both for obj growth and obj ATPM
    #new_model.objective = "ATPM"
    #new_model.reactions.EX_pi_e.lower_bound = -0.14
    #new_model.reactions.EX_pi_e.lower_bound = -10
    new_model.reactions.EX_o2_e.bounds = -5, -5
    new_model.reactions.EX_nh4_e.lower_bound = -1000
    sol = new_model.optimize()
    print(sol)
    print(sol.fluxes['EX_o2_e'])
    print(sol.fluxes['EX_glyc_e'])
    #print('ADK flux', print(sol.fluxes['ADK2']))
    #print(new_model.metabolites.atp_c.summary())
    print(new_model.metabolites.o2_c.summary())

    # limit oxygen consuming complexes
    new_model.reactions.CYTBD.upper_bound = 1
    new_model.reactions.CYTBD2pp.upper_bound = 1
    new_model.reactions.CYTBDpp_1.upper_bound = 1
    new_model.reactions.CYTB_B2.upper_bound = 1
    new_model.reactions.CYO1b.upper_bound = 1
    new_model.reactions.DESAT16.upper_bound = 0
    new_model.reactions.DESAT18.upper_bound = 0
    new_model.reactions.ROO.upper_bound = 0
    new_model.reactions.NOX.upper_bound = 0
    new_model.reactions.GLYCTO1.upper_bound = 0
    new_model.reactions.ASPO1.upper_bound = 0
    new_model.reactions.DALAOX.upper_bound = 0
    new_model.reactions.PYAM5PO.upper_bound = 0
    new_model.reactions.MEHLER.upper_bound = 0
    #after fixing o2 bounds:
    new_model.reactions.FEROpp.upper_bound = 0
    new_model.reactions.GGPTRCO.upper_bound = 0
    #after moving to minimal medium without soy protein
    new_model.reactions.HGNTOR.upper_bound = 0
    new_model.reactions.get_by_id("34HPPOR").upper_bound = 0
    new_model.reactions.get_by_id("FAS181").upper_bound = 0
    new_model.reactions.get_by_id("FAS161").upper_bound = 0
    new_model.reactions.get_by_id("FDMO2").upper_bound = 0
    new_model.reactions.get_by_id("FDMO_1").upper_bound = 0
    new_model.reactions.get_by_id("LYSMO").upper_bound = 0
    new_model.reactions.get_by_id("INOSTO").upper_bound = 0
    new_model.reactions.get_by_id("PYROX").upper_bound = 0
    sol = new_model.optimize()
    print('results after modulating O2-consuming reactions')
    print(sol)
    print(sol.fluxes['EX_o2_e'])
    print(sol.fluxes['EX_glyc_e'])
    print("clav flux is: " , sol.fluxes['EX_clav_e'])
    print(sol.fluxes['EX_pi_e'])
    # print(new_model.metabolites.o2_c.summary())
    # print(new_model.metabolites.o2_p.summary())
    # print(new_model.metabolites.clavam_c.summary())
    # print(new_model.metabolites.clavad_c.summary())
    # print(new_model.metabolites.cmeclav_c.summary())
    # print(new_model.metabolites.fomclav_c.summary())
    # print(new_model.metabolites.hmclav_c.summary())
    # print(new_model.metabolites.fclav_c.summary())
    # print(new_model.metabolites.clavcar_c.summary())
    print("clavam-2-carboxylate flux is: ", sol.fluxes['EX_clavcar_e'])

    #print(new_model.metabolites.atp_c.summary(solution=sol))

    for reac in new_model.reactions:
        if reac.id.startswith('EX_'):
            flux = sol.fluxes[reac.id]
            if flux != 0:
                print(reac.id, flux)
    #print(new_model.metabolites.h_c.summary(solution=sol))
    print(new_model.metabolites.o2_c.summary(solution=sol))
    return None



def make_overview_different_mass_balances(ori_model_file, new_model_file):
    ori_model = cobra.io.read_sbml_model(ori_model_file)
    new_model = cobra.io.read_sbml_model(new_model_file)

    diff_list = []

    for new_reaction in new_model.reactions:
        # find the matching ori reaction
        try:
            ori_reaction = ori_model.reactions.get_by_id(new_reaction.id)
        except KeyError:
            #print(new_reaction.id)
            continue
        ori_meta_dict = ori_reaction.metabolites
        new_meta_dict = new_reaction.metabolites
        #convert the metabolites to the same id's
        upd_ori_meta_dict = {}
        upd_new_meta_dict = {}
        for key in ori_meta_dict:
            first_part_id, second_part_id = key.id.split("[")
            new_id = f"{first_part_id}_{second_part_id[0]}"
            upd_ori_meta_dict[new_id] = ori_meta_dict[key]
        for key in new_meta_dict:
            upd_new_meta_dict[key.id] = new_meta_dict[key]
        if upd_ori_meta_dict != upd_new_meta_dict:
            # only keep the differences between the two dicts
            new_line = [new_reaction.id]
            combined_keys = upd_ori_meta_dict.keys() | upd_new_meta_dict.keys()
            for key in combined_keys:
                try:
                    ori_value = upd_ori_meta_dict[key]
                except KeyError:
                    new_line.append(f"Ori_model lacks {key}:{upd_new_meta_dict[key]}")
                    continue
                try:
                    new_value = upd_new_meta_dict[key]
                except KeyError:
                    new_line.append(f"New_model lacks {key}:{upd_ori_meta_dict[key]}")
                if ori_value != new_value:
                    new_line.append(key)
                    new_line.append(f"ori value: {ori_value}")
                    new_line.append(f"new_value: {new_value}")
            print(new_line)
            diff_list.append(new_line)
        else:
            #print(new_reaction.id)
            pass
    print(len(diff_list))
    return None


def rewrite_gene_ids_ori_model(input_file, blast_result_file, output_files, difference_overview):
    model = cobra.io.read_sbml_model(input_file)
    blast_df = pd.read_csv(blast_result_file, sep="\t")
    list_of_genes = [gene.id for gene in model.genes]
    filtered_blast_df = blast_df[blast_df['qseqid'].isin(list_of_genes)]

    for sorter_column, output_tsv in zip(['pident', 'bitscore'], output_files):

        filtered_blast_df = filtered_blast_df.sort_values(sorter_column, ascending=False)
        geneIDs_matrix = filtered_blast_df.drop_duplicates(subset='qseqid', keep='first')

        geneIDs_matrix.set_index('qseqid', inplace=True)

        geneIDs_matrix_sorted = geneIDs_matrix.sort_index()

        geneIDs_matrix_sorted.to_csv(output_tsv, sep='\t')

    # compare the two data_frames
    data_frame_one = pd.read_csv(output_files[0], sep="\t", index_col=0)
    data_frame_two = pd.read_csv(output_files[1], sep="\t", index_col=0)
    compare_dfs = data_frame_one['sseqid'].compare(data_frame_two['sseqid'])
    print(compare_dfs.shape)
    print(compare_dfs.head())
    compare_dfs.to_csv(difference_overview, sep='\t')

    print(data_frame_one.shape)
    print(data_frame_one['sseqid'].nunique())
    print(data_frame_two.shape)
    print(data_frame_two['sseqid'].nunique())
    return None


def check_gpr_ori_and_new_model(ori_model_file, new_model_file, gene_id_mismatch_file):
    ori_model = cobra.io.read_sbml_model(ori_model_file)
    new_model = cobra.io.read_sbml_model(new_model_file)
    gene_id_matches = pd.read_csv(gene_id_mismatch_file, sep="\t")

    for row in gene_id_matches.iterrows():
        old_gene_id = row[1]['qseqid']
        print(f"starting check for {old_gene_id}")
        wrong_gene_id = row[1]['self']
        correct_gene_id = row[1]['other']
        old_gene_object = ori_model.genes.get_by_id(old_gene_id)
        # print(dir(old_gene_object))
        try:
            new_gene_object = new_model.genes.get_by_id(wrong_gene_id)
            print(f"{wrong_gene_id} is (wrongfully) in the model")
            print(f"reactions for ori model are {old_gene_object.reactions}")
            print(f"reactions for new model are {new_gene_object.reactions}")
            print(f"the correct id should be {correct_gene_id}")
        except KeyError:
            print(f"{wrong_gene_id} is not in the model")
            try:
                new_gene_object = new_model.genes.get_by_id(correct_gene_id)
                print(f"{correct_gene_id} is in the model instead")
                print(f"reactions for ori model are {old_gene_object.reactions}")
                print(f"reactions for new model are {new_gene_object.reactions}")
            except KeyError:
                print(f"{correct_gene_id} is also not in the model")
    return None


def quick_checks_updated_model(model_file_upd, old_model_file):
    upd_model = cobra.io.read_sbml_model(model_file_upd)
    old_model = cobra.io.read_sbml_model(old_model_file)

    # for reac in ["NADH16","NADH2r","NADH4","NADH5","NADH6"]:
    #     print(old_model.reactions.get_by_id(reac).gene_reaction_rule)
    #     print(model.reactions.get_by_id(reac).gene_reaction_rule)

    # for gene in model.genes:
    #     if not gene.id.startswith("SCLAV"):
    #         print(gene.id)

    print("old medium")
    print(old_model.medium)
    print("new medium")
    print(upd_model.medium)

    upd_model.reactions.EX_pi_e.lower_bound = -0.14
    old_model.reactions.EX_pi_e.lower_bound = -0.14

    old_sol = old_model.optimize()
    upd_sol = upd_model.optimize()

    print(old_sol.fluxes["EX_nh4_e"])
    print(upd_sol.fluxes["EX_nh4_e"])

    for reac in upd_model.reactions:
        if reac.id.startswith('EX_'):
            flux = upd_sol.fluxes[reac.id]
            if flux != 0:
                print(reac.id, flux)
    return None


def remove_soy_reaction(in_model_file, out_model_name):
    model = cobra.io.read_sbml_model(in_model_file)
    model.reactions.SOYPEPT.remove_from_model()
    cobra.io.write_sbml_model(model, out_model_name)
    return None

if __name__ == '__main__':
    ori_mod_file = "Model_files/iDG1237.xml"
    first_mod_file = "Model_files/NewModel_CleanedFormula_Hbalance3.xml"
    #new_mod_file = "Model_files/iDG1237_updated.xml"
    new_mod_file = "Model_files/iDG1237_updated_with_uniprot.xml"
    upd_mod_file = "Model_files/iDG1237_updated_with_uniprot_v2.xml"
    #compare_models(ori_mod_file, new_mod_file)
    compare_models(ori_mod_file, upd_mod_file)

    #try_out_clav_production_cleaned_formula_model(first_mod_file)
    #try_out_clav_production_updated_uniprot_model(new_mod_file)

    #make_overview_different_mass_balances(ori_mod_file, new_mod_file)
    blast_result = "input_files/blast_sclav_results.tsv"
    gene_id_file_pident = "input_files/new_gene_ids_lookup_table_pident.tsv"
    gene_id_file_bitscore = "input_files/new_gene_ids_lookup_table_bitscore.tsv"
    gene_id_diffs_file = "input_files/gene_id_mismatches_pident_to_bitscore.tsv"
    # rewrite_gene_ids_ori_model(ori_mod_file, blast_result, [gene_id_file_pident, gene_id_file_bitscore],
    #                            gene_id_diffs_file)

    # check_gpr_ori_and_new_model(ori_mod_file, new_mod_file, gene_id_diffs_file)

    #quick_checks_updated_model(new_mod_file, ori_mod_file)
    # quick_checks_updated_model(new_mod_file, first_mod_file)

    #remove_soy_reaction(new_mod_file, upd_mod_file)

