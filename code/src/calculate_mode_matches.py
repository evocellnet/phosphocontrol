"""
"""

from sys import argv
from pathlib import Path
import subprocess
from multiprocessing import Pool
import re

from itertools import combinations
import numpy as np
from scipy.stats import spearmanr, pearsonr
from scipy.optimize import linear_sum_assignment
from sklearn.metrics.pairwise import cosine_similarity
import pandas as pd
from tqdm import tqdm



def call_r_script(modes_file_arg, out_arg):
    """
    Call R script to dump the eigenvectors
    
    Arguments
    ---------
    modes_file_arg: path to NMA Rdata file
    out_arg: path for output csvs

    Returns none

    """

    call = ['Rscript','write_eigenvector_matrix.R',modes_file_arg, out_arg]

    try:
        subprocess.check_call(call)
    except subprocess.CalledProcessError as error:
        raise ValueError(f'Rscript failed: {error}')
    except OSError:
        raise OSError("Rscript executable not found in PATH")



def define_prot_to_status(annotation_df):
    """
    """
    prot_to_status = {}
    for _, row in annotation_df.iterrows():
        prot_id = row["Full ID"]
        status = row["Status"]
        prot_to_status[prot_id] = status
    return prot_to_status

def compute_pairwise_similarities(eigenvector_mtx_a, eigenvector_mtx_b):
    """
    Calculate the cosine similarities between eigenvectors

    Arguments
    ---------
    eigenvector_mtx_X: matrix where the columns are eigenvectors

    Returns
    -------
    similarity_matrix: MxM matrix of cosine similarities

    """
    similarity_matrix = cosine_similarity(eigenvector_mtx_a.T, eigenvector_mtx_b.T)
    return similarity_matrix

def solve_assignment(cost_matrix):
    """
    Find best match in the cost matrix using the Hungarian algorithm

    Arguments
    ---------
    cost_matrix: array-like

    Returns
    -------
    row_ind: array of row indices
    col_ind: array of column indices
    cost: float, total cost of the assignment

    row_ind and col_ind provides the optimal assignment. E.g.
    array([0, 1, 2]), array([1, 0, 2]))
    row 0 matches column 1
    row 1 matches column 0
    row 2 matches column 2
    """
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    return row_ind, col_ind, cost_matrix[row_ind, col_ind].sum()

def get_matching_similarities(similarity_matrix, row_ind, col_ind):
    """
    Retrieve similarities between matching modes
    """
    cosine_similarities = []
    for idx, i in enumerate(row_ind):
        j = col_ind[idx]
        cosine_sim = similarity_matrix[i,j]
        cosine_similarities.append(cosine_sim)
    return cosine_similarities


def extract_pdb_code(filename):
    """
    Function to extract PDB code from the eigenvector filenames
    """
    match = re.search(r"eigenvectors_(.*?)\.csv", filename)
    if match:
        return match.group(1)
    return None


def process_valid_subdir(args):
    """
    """

    valid_subdir, annot_path, input_path, out_path = args
    psite = valid_subdir.stem

    psite_out_path = out_path / psite
    psite_out_path.mkdir(exist_ok=True)

    matrices_path = psite_out_path / "eigenvector_matrices"
    matrices_path.mkdir(exist_ok=True)

    # sims_path = psite_out_path / "similarity_matrices"
    # sims_path.mkdir(exist_ok=True)

    annot_df = pd.read_csv(annot_path / f"{psite}.csv")
    prot_to_status = define_prot_to_status(annot_df)
    
    # Dump eigenvector matrices
    mode_file = input_path / psite / "modes.RData.gz"
    call_r_script(mode_file, matrices_path)
    eigenvector_files_list = list(matrices_path.glob("*.csv"))
    eigenvector_combs = list(combinations(eigenvector_files_list, 2))

    # Calculate all pairwise similarities between normal modes
    within_p_similarities = []
    within_np_similarities = []
    between_similarities = []
    corr_df = []

    for comb in eigenvector_combs:
        matrix_a_f, matrix_b_f = comb
        
        pdb_code_a = extract_pdb_code(str(matrix_a_f.parts[-1]))
        pdb_code_b = extract_pdb_code(str(matrix_b_f.parts[-1]))
        status_a = prot_to_status[pdb_code_a]
        status_b = prot_to_status[pdb_code_b]

        matrix_a = pd.read_csv(matrix_a_f)
        matrix_b = pd.read_csv(matrix_b_f)
        similarity_matrix = compute_pairwise_similarities(matrix_a, matrix_b)
        
        # # For debugging
        # sim_mtx_df = pd.DataFrame(similarity_matrix)
        # sim_mtx_df.to_csv(sims_path / f"matrix_{pdb_code_a}_{pdb_code_b}.csv")

        # Compute the optimal matching between normal modes
        cost_matrix = 1 - similarity_matrix
        row_ind, col_ind, cost = solve_assignment(cost_matrix)
        best_matching_similarities = get_matching_similarities(similarity_matrix, row_ind, col_ind)
        spearman_corr, spearman_pval = spearmanr(row_ind, col_ind)
        pearson_corr, pearson_pval = pearsonr(row_ind, col_ind)

        if status_a == "Phosphorylated" and status_b == "Phosphorylated":
            #print(psite, pdb_code_a,status_a, pdb_code_b, status_b, "within p")
            within_p_similarities.append(best_matching_similarities)
            corr_df.append([spearman_corr, spearman_pval, pearson_corr, pearson_pval, 'within_phospho'])
        elif status_a == "Non-phosphorylated" and status_b == "Non-phosphorylated":
            #print(psite, pdb_code_a,status_a, pdb_code_b,status_b, "within np")
            within_np_similarities.append(best_matching_similarities)
            corr_df.append([spearman_corr, spearman_pval, pearson_corr, pearson_pval, 'within_nonphospho'])
        else:
            #print(psite, pdb_code_a,status_a, pdb_code_b,status_b, "between")
            between_similarities.append(best_matching_similarities)
            corr_df.append([spearman_corr, spearman_pval, pearson_corr, pearson_pval, 'between_groups'])
    
    columns = [f"mode_{str(i)}" for i in range(0,len(best_matching_similarities))]

    within_p_similarities = pd.DataFrame(within_p_similarities,columns=columns)
    within_np_similarities = pd.DataFrame(within_np_similarities,columns=columns)
    between_similarities = pd.DataFrame(between_similarities,columns=columns)

    within_p_similarities.to_csv(psite_out_path / "within_phospho_modes.csv")
    within_np_similarities.to_csv(psite_out_path / "within_nonphospho_modes.csv")
    between_similarities.to_csv(psite_out_path / "between_modes.csv")

    corr_df = pd.DataFrame(corr_df, columns=['spearman_rho','spearman_pval',
                                             'pearsonr','pearson_pval','status'])
    corr_df.to_csv(psite_out_path / "corr_df.csv")


def main(valid_subdirs, annot_path, input_path, out_path, num_workers=4):
    """
    Parallelizes process_valid_subdir 

    Arguments
    ---------
    valid_subdirs: list of valid directories to take input from
    annot_path: directory containing chain annotation dataframes
    input_path: base input directory
    out_path: base output directory
    
    Returns
    -------
    None, writes to disk
    """

    # Prepare arguments for parallel processing
    args = [(valid_subdir, annot_path, input_path, out_path) for valid_subdir in valid_subdirs]

    with Pool(num_workers) as pool:
        list(tqdm(pool.imap(process_valid_subdir, args), total=len(valid_subdirs)))
    pool.join()


if __name__ == "__main__":

    #########
    # Setup #
    #########

    num_workers = int(argv[1])

    out_path = Path("../../results/new_nma_mode_matching")
    out_path.mkdir(exist_ok=True)
    input_path = Path("../../results/new_nma_analysis_per_psite")
    #input_path = Path("/cluster/work/beltrao/mcorrea/new_nma_analysis_per_psite")
    annot_path = Path("../../data/processed/pdb_pairs/chains_by_protein/annotation_per_psite")
    subdirs = [x for x in input_path.iterdir() if x.is_dir()]

    # Get phosphosites for which we will do the analysis
    flucts_df_nodups_f = Path("../../notebooks/global_dynamics_overview/fluctuations_df_nodups.csv")
    flucts_df_nodups = pd.read_csv(flucts_df_nodups_f, index_col=0)
    valid_psites = list(flucts_df_nodups["phosphosite"])

    # Get examples for which we have the NMA saved
    # (i.e. those where it could be done)
    valid_subdirs = []
    mode_files = []
    for subdir in subdirs:
        psite = subdir.stem
        mode_file = subdir / "modes.RData.gz"
        if mode_file.is_file() and psite in valid_psites:
            valid_subdirs.append(subdir)
            mode_files.append(mode_file)

    ###################
    # Match NMA modes #
    ###################

    main(valid_subdirs[:20], annot_path, input_path, out_path, num_workers=num_workers)
    
    print("Done!")


