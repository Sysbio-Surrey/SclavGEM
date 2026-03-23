import cobra
import pandas as pd
import ast


def delete_metabolites(model, metabolite_file):
    metabolite_df = pd.read_csv(metabolite_file, sep='\t', header=None, names=['meta_id'])
    for meta_id in metabolite_df['meta_id']:
        model.metabolites.get_by_id(meta_id).remove_from_model()
    return model


def rework_metabolite_ids(model):
    """metabolite ids are reworked from the metabolite[comp] format to metabolite_comp

    :param model:
    :return:
    """
    for metabolite in model.metabolites:
        old_id_part_one, old_id_part_two = metabolite.id.split("[")
        new_id = f"{old_id_part_one}_{old_id_part_two[0]}"
        metabolite.id = new_id
    return model


def update_metabolite_formulas(model, updates_file):
    formula_df = pd.read_csv(updates_file)
    changed_metabolites = formula_df["metabolite_id"].values
    for metabolite in model.metabolites:
        if metabolite.id in changed_metabolites:
            new_formula = formula_df.loc[formula_df["metabolite_id"] == metabolite.id, "new_formula"].iloc[0]
            metabolite.formula = new_formula
    return model


def add_new_metabolites(model, add_meta_file):
    meta_df = pd.read_csv(add_meta_file, sep="\t")
    for idx, row in meta_df.iterrows():
        new_meta = cobra.Metabolite(
            id=row['id'],
            name=row['name'],
            formula=row['formula'],
            charge=row['charge'],
            compartment=row['compartment'],
        )
        new_meta.annotation = ast.literal_eval(row['annotation'])
        new_meta.notes = ast.literal_eval(row['notes'])

        model.add_metabolites([new_meta])
    return model


if __name__ == "__main__":
    # this script summarizes all the changes to the model metabolites
    input_model = "model/iDG1237.xml"
    # input files
    deleted_meta_file = "data/data_to_update_published_model/Metabolites_deleted_from_ori_model.txt"
    update_metabolite_formula_file = "data/data_to_update_published_model/changed_metabolite_formulas.csv"
    new_meta_file = "data/data_to_update_published_model/new_metabolites.tsv"

    # read model
    input_model_obj = cobra.io.read_sbml_model(input_model)

    # delete metabolites
    delete_metabolites(input_model_obj, deleted_meta_file)
    # update model id's
    upd_model = rework_metabolite_ids(input_model_obj)
    # update metabolite formula's
    update_metabolite_formulas(upd_model, update_metabolite_formula_file)
    # add new metabolites
    add_new_metabolites(upd_model, new_meta_file)

    #missing functions: add metabolites
    #remove metabolites...


