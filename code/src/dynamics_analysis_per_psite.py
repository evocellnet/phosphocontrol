"""
Call R NMA modelling script for all proteins
The analysis is done per phosphosite
"""
import pandas as pd
from pathlib import Path
import subprocess
from tqdm import tqdm


def prepare_annotation_dfs(df, out_path, structures_path):

    psite_to_annot_df = {}
    columns = ["Protein ID","PDB ID","Full ID","Path","Status","color"]
    unique_psites = list(df["PHOSPHOSITE"].unique())

    for unique_psite in unique_psites:

        annot_df = []

        uniprot_id, res_idx = unique_psite.split("_")
        subset_df = df.loc[df["PHOSPHOSITE"]==unique_psite]
        nonphospho_names = list(subset_df["RENAMED_ENTITY_ONE"].unique())
        phospho_names = list(subset_df["RENAMED_ENTITY_TWO"].unique())

        nonphospho_rows = produce_rows(nonphospho_names, uniprot_id, structures_path, False)
        phospho_rows = produce_rows(phospho_names, uniprot_id, structures_path, True)

        annot_df = nonphospho_rows + phospho_rows
        annot_df = pd.DataFrame(annot_df, columns = columns)
        annot_df_path = out_path / f"{unique_psite}.csv"
        annot_df.to_csv(annot_df_path ,index=False)
        psite_to_annot_df[unique_psite] = annot_df_path
    
    return psite_to_annot_df


def produce_rows(names, uniprot_id, structures_path, is_phosphorylated):
    rows = []
    for name in names:
        pdb_id, _ = name.split("_")
        chain_path = structures_path / uniprot_id / f"{name}.pdb"
        if is_phosphorylated:
            color = "orange"
            status = "Phosphorylated"
        else:
            color = "blue"
            status = "Non-phosphorylated"
        
        row = [uniprot_id, pdb_id, name, chain_path, status, color]
        rows.append(row)
    return rows

def call_r_script(input_arg, annot_arg, out_arg , stdout_f):
    """
    """
    call = ['Rscript','nma_analysis.R',input_arg, annot_arg, out_arg]

    #print(call)
    with open(stdout_f, 'w') as f:
        try:
            subprocess.check_call(call, stdout=f)
        except subprocess.CalledProcessError as error:
            raise ValueError(f'Rscript failed: {error}')
        except OSError:
            raise OSError("Rscript executable not found in PATH")


if __name__ == '__main__':

    #########
    # Setup #

    data_path = Path("../../data/processed/pdb_pairs")
    df = pd.read_csv(data_path / "filtered_df.csv")

    chains_path = Path(data_path / "chains_by_protein")
    annotation_path = chains_path / "annotation_per_psite"
    annotation_path.mkdir(exist_ok=True)
    base_output = Path("../../results/nma_analysis_per_psite")
    #base_output = Path("/cluster/work/beltrao/mcorrea/new_nma_analysis_per_psite")
    base_output.mkdir(exist_ok=True)
    # Create annotation dataframes
    psite_to_annot_df = prepare_annotation_dfs(df, annotation_path, chains_path)

    ###############
    # Call script #
    ###############

    for psite, annot_df_path in tqdm(psite_to_annot_df.items()):
        print(f"Now on {psite}")
        uniprot_id, _ = psite.split("_")
        subdir = chains_path / uniprot_id
        output_path = base_output / psite
        stdout_file = base_output / f"{psite}.txt"

        try:
            call_r_script(str(subdir), str(annot_df_path), str(output_path), str(stdout_file))
        except ValueError:
            print(f"Error on {psite}")
            continue
    
    print("Et voila!")

