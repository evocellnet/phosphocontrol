"""
Script to mine AUTHOR field data and wrangle it
"""

from pathlib import Path
import pandas as pd
from collections import Counter
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from tqdm import tqdm


def get_structure_to_lead_author(valid_cif_files):
    """
    Create dictionary linking a structure to its lead author
    """

    structure_to_lead_author = {}
    for cif_file in tqdm(valid_cif_files):
        mmcif_dict = MMCIF2Dict(cif_file)
        try:
            authors = mmcif_dict["_citation_author.name"]
            lead_author = authors[-1]
            structure_to_lead_author[cif_file.stem.lower()] = lead_author
        except KeyError:
            print(f"File {cif_file} missing author name; adding placeholders")
            structure_to_lead_author[cif_file.stem.lower()] = "ANONYMOUS"

    # Hack: manually add author name https://www.rcsb.org/structure/3UVR
    structure_to_lead_author["3uvr"] = "Rauh, D."

    return structure_to_lead_author


def prepare_model_df(subdirs):
    """
    Create a dataframe that can be used as input for a linear model.

    Returns
    -------
    glm_df: pandas DataFrame
            Contains the following columns:
            structure_pair: PDB ID of the two structures, joined by $
            rmsd: RMSD between the two structures
            is_phospho_a, is_phospho_b: whether the structures are phosphorylated
            same_lead_author: whether the two structures share the same lead author
    """

    glm_df = []
    print('Preparing the glm dataframe!')
    columns = ["structure_pair","rmsd","is_phospho_a","is_phospho_b","same_lead_author"]
    for subdir in tqdm(subdirs):
        print(f"Now on {subdir}")
        rmsd_df = pd.read_csv(subdir / "rmsds_df.csv",index_col=0,
                              dtype={'PDB_ID_A':'string',
                                     'PDB_ID_B': 'string'})
        
        rmsd_df["ID_A"] = rmsd_df["PDB_ID_A"] + "_" + rmsd_df["Chain_A"]
        rmsd_df["ID_B"] = rmsd_df["PDB_ID_B"] + "_" + rmsd_df["Chain_B"]
        
        print(f'Iterating over RMSD dataframe for {subdir}')
        for _, row in tqdm(rmsd_df.iterrows()):
            chain_a = row["ID_A"]
            chain_b = row["ID_B"]
            pdb_code_a = chain_a.split("_")[0]
            pdb_code_b = chain_b.split("_")[0]

            structure_pair = "$".join([chain_a, chain_b])
            rmsd = row["RMSD"]
            group = row["Group"]
            
            lead_author_a = structure_to_lead_author['lead_author'][pdb_code_a]
            lead_author_b = structure_to_lead_author['lead_author'][pdb_code_b]
            
            if lead_author_a == lead_author_b:
                same_lead_author = 1
            else:
                same_lead_author = 0
            
            if group == "between_groups":
                is_phospho_a = 0
                is_phospho_b = 1
            elif group == "within_phospho":
                is_phospho_a = 1
                is_phospho_b = 1
            elif group == "within_nonphospho":
                is_phospho_a = 0
                is_phospho_b = 0
            else:
                raise ValueError(f"Invalid group value {group}")
            row = [structure_pair, rmsd, is_phospho_a, is_phospho_b, same_lead_author]
            glm_df.append(row)
    
    glm_df = pd.DataFrame(glm_df,columns=columns)
    return glm_df


if __name__ == "__main__":

    #########
    # Setup #
    #########

    raw_pdbs = Path("../../data/raw/pdb_pairs/pdb_files/pdb")
    pdb_files = list(raw_pdbs.glob("*.pdb"))

    raw_cifs = Path("../../data/raw/pdb_pairs/pdb_files/cif")
    cif_files = list(raw_cifs.glob("*.cif"))

    # Get the structures that we actually use
    df = pd.read_csv("../../data/processed/pdb_pairs/filtered_df.csv")
    pdbs_a = set(df["ENTRY_ID_ONE"])
    pdbs_b = set(df["ENTRY_ID_TWO"])
    valid_pdbs = pdbs_a.union(pdbs_b)
    print(f"Valid number of PDB IDs: {len(valid_pdbs)}")

    valid_cif_files = [x for x in cif_files if x.stem.lower() in valid_pdbs]
    print(f"{len(valid_cif_files)} valid CIF files")

    ###################################
    # Associate structures to authors #
    ###################################

    print('Associating structures to authors...')
    structure_to_lead_author = get_structure_to_lead_author(valid_cif_files)
    author_counter = Counter(structure_to_lead_author.values())
    author_df = pd.DataFrame.from_dict(structure_to_lead_author,orient='index',
                                       columns=['lead_author'])
    author_df.to_csv("../../data/processed/pdb_pairs/author_df.csv")

    # Just upload the file and turn it into a dict
    # Otherwise we need to upload 19 gb to the cluster
    #author_df_path = Path("../../data/processed/pdb_pairs/author_df.csv")
    #author_df = pd.read_csv(author_df_path, index_col=0)
    #structure_to_lead_author = author_df.to_dict()

    #############################################
    # Wrangle data for linear model #
    #############################################

    rmsds_path = Path("../../results/rmsds")
    subdirs = [x for x in rmsds_path.iterdir() if x.is_dir()]

    glm_df = prepare_model_df(subdirs)
    print('Preparing dataframe...')
    glm_df.to_csv("../../data/processed/glm_control_input/glm_df.csv")

    phospho_df = glm_df.loc[(glm_df["is_phospho_a"]==1) & (glm_df["is_phospho_b"]==1)]
    phospho_df.to_csv("../../data/processed/glm_control_input/phospho_df.csv")
    nonphospho_df = glm_df.loc[(glm_df["is_phospho_a"]==0) & (glm_df["is_phospho_b"]==0)]
    nonphospho_df.to_csv("../../data/processed/glm_control_input/nonphospho_df.csv")

    print('Done!')