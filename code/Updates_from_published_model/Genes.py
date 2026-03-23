import cobra
import pandas as pd


def update_gene_ids(model, gene_id_table):
    gene_id_lookup_df = pd.read_csv(gene_id_table, sep="\t")
    for gene in model.genes:
        try:
            new_gene_id = gene_id_lookup_df.loc[gene_id_lookup_df["qseqid"] == gene.id, "sseqid"].iloc[0]
            gene.id = new_gene_id
        except IndexError:
            print(gene.id)
            continue
    return model


def delete_genes(model, delete_list_csv):


if __name__ == "__main__":
    model_file = "../Model_files/iDG1237.xml"
    model = cobra.io.read_sbml_model(model_file)
    #step 1: update gene-ids
    gene_id_lookup_table = "../input_files/new_gene_ids_lookup_table_bitscore.tsv"
    update_gene_ids(model, gene_id_lookup_table)
    #step 2: remove genes

    #step 3: update gpr-rules