"""
Quick adaptation of all_global_rmsds to better deal with cases with massive
amounts of structure pairs
"""

import argparse
from pathlib import Path
from itertools import combinations, product
from multiprocessing import Pool
from warnings import warn

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from tqdm import tqdm

from rmsd import rmsd
from all_global_rmsds import get_unique_phosphosites, get_domain_phosphosites, get_domain_coverage

global CHUNK_SIZE
CHUNK_SIZE = int(5e3)

def digest_args():
    """
    Parse script arguments
    """

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input_pairs", type=Path,
                    help="Structure pair dataframe")
    ap.add_argument("--structures", type=Path, help="Path to structural data")
    ap.add_argument("--domain_mode", type=bool,
                    help="Whether to use only the domains where phosphosites\
                          are, rather than whole chains")
    ap.add_argument("--domains", type=Path,
                    help="Phosphosite domain location data", required=False)
    ap.add_argument("--domain_cov",type=Path,
                    help='Path to domain quality information')
    ap.add_argument("--output", type=Path, help="Output path")
    ap.add_argument("--n_cpus", type=int, default=4,
                    help="Number of parallel processes")

    args = ap.parse_args()
    return args


def compare_pair(combination, gap_open=-10, gap_extend=-0.5):
    """
    Calculate RMSDs for a given combination of structures
    """
    ref, mob = combination

    ref_pdb, ref_chain_id, ref_res, reference_path = ref
    mob_pdb, mob_chain_id, mob_res, mobile_path = mob

    parser = PDBParser(QUIET=True)
    ref_mod = parser.get_structure("X", str(reference_path))
    mob_mod = parser.get_structure("Y", str(mobile_path))

    # Deal with NMR ensembles
    # Make all possible combinations of conformers in the files and compute
    # RMSDs on them
    ref_models = [ref_mod[idx][ref_chain_id] for idx, _ in enumerate(ref_mod)]
    mob_models = [mob_mod[idx][mob_chain_id] for idx, _ in enumerate(mob_mod)]
    model_combs = product(ref_models, mob_models)

    model_rmsds = []
    for model_comb in model_combs:
        reference, mobile = model_comb
        superimposition_rmsd = rmsd(reference, mobile, gap_open, gap_extend)
        model_rmsds.append(superimposition_rmsd)

    # If we have an NMR model, the final RMSD between the two files will be
    # the median of all RMSDs between conformers
    median_rmsd = np.median(model_rmsds)

    row = [ref_pdb, ref_chain_id, ref_res,
           mob_pdb, mob_chain_id, mob_res,
           median_rmsd]

    return row


def parallel_rmsds(combinations, group, num_workers):
    """
    """
    with Pool(num_workers) as pool:
        rows = pool.map(
            compare_pair, tqdm(combinations, total=len(combinations)))
    pool.join()

    header = ["PDB_ID_A", "Chain_A", "Residue_A",
              "PDB_ID_B", "Chain_B", "Residue_B",
              "RMSD", "Group"]
    rmsd_df = []
    for row in rows:
        row.append(group)
        rmsd_df.append(row)
    rmsd_df = pd.DataFrame(rmsd_df, columns=header)
    return rmsd_df


def chunk_list(lst, chunk_size=CHUNK_SIZE):
    """Divide a list into chunks of specified size."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

if __name__ == "__main__":

    #########
    # Setup #
    #########

    args = digest_args()
    df = pd.read_csv(args.input_pairs)
    chains_path = args.structures
    num_workers = args.n_cpus
    domain_mode = args.domain_mode
    output_path = args.output
    output_path.mkdir(exist_ok=True)

    # To use the extracted domains, we need to use the file that links each
    # phosphosite to a particular domain
    domain_mode = args.domain_mode
    if domain_mode:
        psite_to_domain = pd.read_csv(args.domains)
        domain_psites, psite_to_pfam_id = get_domain_phosphosites(
            psite_to_domain)
        domain_coverage_df = pd.read_csv(args.domain_cov)
        domain_coverage = get_domain_coverage(domain_coverage_df)

    ####################################################################
    # Organize things according to unique combinations of phosphosites #
    # and unique combinations of:                                      #
    # * phosphorylated vs non-phosphorylated structures                #
    # * phosphorylated vs phosphorylated structures                    #
    # * non-phosphorylated vs non-phosphorylated structures            #
    ####################################################################

    if domain_mode:
        # If we are using only the domains, rather than full chains, do not take
        # into account phosphosites outside detected Pfam domains
        psite_dict = get_unique_phosphosites(
            df, chains_path, domain_psites, psite_to_pfam_id,domain_coverage)
    else:
        psite_dict = get_unique_phosphosites(df, chains_path)


    unique_phosphosites = list(psite_dict.keys())
    print(f"Found {len(unique_phosphosites)} unique phosphosites")
    phospho_combs = {}
    nonphospho_combs = {}

    print("Finding possible pairs of structures...")
    for k, v in tqdm(psite_dict.items()):
        phospho_struct = list(set([x[1] for x in v]))
        nonphospho_struct = list(set([x[0] for x in v]))

        phospho_combs[k] = list(combinations(phospho_struct, 2))
        nonphospho_combs[k] = list(combinations(nonphospho_struct, 2))

    ############################
    # Compute all global RMSDs #
    ############################

    input_list = []
    print("Preparing input...")
    for phosphosite in psite_dict:
        between_combs = psite_dict[phosphosite]
        p_combs = phospho_combs[phosphosite]
        np_combs = nonphospho_combs[phosphosite]
        input_list.append(
            (phosphosite, between_combs, p_combs, np_combs, output_path))
        
    for item in tqdm(input_list):
        psite, between_combs, p_combs, np_combs, output_path = item
        print(f"Now on {psite}...")
        
        psite_path = output_path / psite
        psite_path.mkdir(exist_ok=True)
        rmsd_df_path = psite_path / "rmsds_df.csv"
        # Avoid re-running examples
        is_already_done = rmsd_df_path.is_file()
        if not is_already_done:
            print(f"Now on phosphosite {psite}...")

            # Divide combinations into chunks in order to not run out of memory
            # Write dataframes at each step, then join them at the end
            print("Computing between group RMSDs...")
            i = 0
            for chunk in chunk_list(between_combs):
                chunk_file = psite_path / f"between_combs_{i}.csv"
                if not chunk_file.is_file():
                    between_combs_df = parallel_rmsds(chunk, "between_groups", num_workers)
                    between_combs_df.to_csv(chunk_file)
                    i += 1

            print("Computing within phosphorylated structures RMSDs...")
            j = 0
            for chunk in chunk_list(p_combs):
                chunk_file = psite_path / f"within_phospho_{j}.csv"
                if not chunk_file.is_file():
                    within_phospho_df = parallel_rmsds(chunk, "within_phospho", num_workers)
                    within_phospho_df.to_csv(chunk_file)
                    j += 1

            print("Computing within non-phoshorylated structures RMSDs...")
            k = 0
            for chunk in chunk_list(np_combs):
                chunk_file = psite_path / f"within_nonphospho_{k}.csv"
                if not chunk_file.is_file():
                    within_nonphospho_df = parallel_rmsds(chunk, "within_nonphospho", num_workers)
                    within_nonphospho_df.to_csv(chunk_file)
                    k += 1

        else:
            warn('Already done; skipping {psite}...',RuntimeWarning)
            continue

        csv_files = [f for f in psite_path.glob("*.csv")]
        dfs = []
        for csv_file in csv_files:
            df = pd.read_csv(csv_file)
            dfs.append(df)

        # Concatenate all the dataframes into a single dataframe
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df.to_csv(rmsd_df_path, index=False)


    print("Finished!")


