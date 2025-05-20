"""
Script to calculate all vs all global RMSDs per phosphosite
For each phosphosite, they are calculated by:
 - All phosphorylated versus non-phosphorylated
 - All phosphorylated versus non-phosphorylated
 - All non-phosphorylated versus non-phosphorylated

This forms the basis of the analyses shown in Fig. 1 concerning RMSD.

Usage:
python all_global_rmsds.py 
    --input_pairs
    --structures
    --domain_mode TRUE/FALSE
    --domains 
    --domain_cov
    --output OUTPUT_PATH
    --n_cpus INTEGER
"""

import os
import matplotlib as mpl
mpl.use('Agg') # Allows running script on cluster
import argparse
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from itertools import combinations

import pandas as pd
import numpy as np
from tqdm import tqdm
from itertools import product

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'

from Bio.PDB import PDBParser
from rmsd import rmsd
from warnings import warn


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
    # TODO: should probably call this something else like "domain_quality"
    ap.add_argument("--domain_cov",type=Path,
                    help='Path to domain quality information')
    ap.add_argument("--output", type=Path, help="Output path")
    ap.add_argument("--n_cpus", type=int, default=4,
                    help="Number of parallel processes")

    args = ap.parse_args()
    return args


def is_file_empty(file_path):
    return (os.stat(str(file_path)).st_size == 0)


def get_unique_phosphosites(df, chains_path, domain_psites=None,
                            psite_to_pfam_id=None, domain_coverage=None,
                            min_cov_threshold=.75,to_skip=[]):
    """
    Collect information per phosphosite.
    A unique phosphosite is defined by Uniprot ID_residue index

    Returns
    -------
    psite_dict: defaultdict of lists, {phosphosite: [(non-phospho tuple, phospho tuple)]}
                Each tuple consists of the PDB ID, chain ID, residue index and file path

    """
    psite_dict = defaultdict(list)
    for idx, row in df.iterrows():
        uniprot_id = row["UNIPROT"]
        if uniprot_id in to_skip:
            continue
        res_idx = row["UNIPROT_RESIDUE_NUMBER"]
        psite = "_".join([uniprot_id, str(res_idx)])

        pdb_a = row["ENTRY_ID_ONE"]
        chain_a = row["RENAMED_ASYM_ID_ONE"]
        res_a = row["AUTH_SEQ_ID_ONE"]

        pdb_b = row["ENTRY_ID_TWO"]
        chain_b = row["RENAMED_ASYM_ID_TWO"]
        res_b = row["AUTH_SEQ_ID_TWO"]

        if domain_psites is None:
            file_a = chains_path / f"{pdb_a}_{chain_a}.pdb"
            file_b = chains_path / f"{pdb_b}_{chain_b}.pdb"
        else:
            try:
                domain_id = psite_to_pfam_id[psite]
                id_prot_a = f"{pdb_a}_{chain_a}_{domain_id}"
                id_prot_b = f"{pdb_b}_{chain_b}_{domain_id}"
                coverage_a = domain_coverage[id_prot_a]
                coverage_b = domain_coverage[id_prot_b]
                file_a = chains_path / f"{id_prot_a}.pdb"
                file_b = chains_path / f"{id_prot_b}.pdb"
                # Skip comparison if coverage is too low
                if coverage_a < min_cov_threshold or coverage_b < min_cov_threshold:
                    continue
            except KeyError:
                # Skip phosphosites not in Pfam domains
                continue

        if (file_a.exists() and file_b.exists()) and not (is_file_empty(file_a) or is_file_empty(file_b)):
            tup_a = (pdb_a, chain_a, res_a, file_a)
            tup_b = (pdb_b, chain_b, res_b, file_b)
            psite_dict[psite].append((tup_a, tup_b))
        else:
            # Safeguard against empty/absent files
            # A domain may not actually be present in a structure; if we are
            # comparing only domains, this will lead to an empty or non-existent
            # file. We need to ignore those
            warn(f"Empty or absent files: {file_a},{file_b}",RuntimeWarning)
            continue

    return psite_dict


def get_domain_phosphosites(psite_to_domain):
    """
    Get the set of phosphosites that are located within Pfam domains, and the
    ID of the domain associated to the phosphosite

    Returns
    -------
    domain_psites: set, phosphosites found in the dataset (Uniprot ID_residue index)
    psite_to_pfam_id: dict, {psite: Pfam domain ID}
    """

    domain_psites = set()
    psite_to_pfam_id = {}

    for _, row in psite_to_domain.iterrows():
        psite = row["Phosphosite"]
        domain_id = row["Pfam ID"]
        domain_psites.add(psite)
        psite_to_pfam_id[psite] = domain_id

    return domain_psites, psite_to_pfam_id

def get_domain_coverage(df):
    """
    """
    df['unique_id'] = df['PDB ID'] + "_" + df['Chain ID'] + "_" + df['Pfam ID']
    unique_ids = list(df['unique_id'].astype(str))
    coverages = list(df['Coverage'])
    domain_coverage = {unique_id:coverages[idx] for idx, unique_id in enumerate(unique_ids)}

    return domain_coverage

def compute_all_rmsds(psite, between_combs, p_combs, dp_combs, output_path):
    """
    Compute all comparisons per phosphosite
    """
    psite_path = output_path / psite
    psite_path.mkdir(exist_ok=True)
    rmsd_df_path = psite_path / "rmsds_df.csv"

    if not rmsd_df_path.is_file():

        print(f"Now computing RMSDs for {psite}...")
        print("Combinations between groups...")
        between_groups = run_comparison(
            between_combs, "between_groups")
        print("Combinations within phospho...")
        within_phospho = run_comparison(
            p_combs, "within_phospho")
        print("Combinations within non-phospho...")
        within_nonphospho = run_comparison(
            dp_combs, "within_nonphospho")
            # Write output
        header = ["PDB_ID_A", "Chain_A", "Residue_A",
              "PDB_ID_B", "Chain_B", "Residue_B",
              "RMSD", "Group"]

        rmsd_df = between_groups + within_nonphospho + within_phospho
        rmsd_df = pd.DataFrame(rmsd_df, columns=header)
        rmsd_df.to_csv(str(rmsd_df_path))
        plot_histogram(rmsd_df, psite_path)
        plot_boxplot(rmsd_df, psite_path)
    else:
        print("Already calculated, skipping...")



def run_comparison(combinations, group, gap_open=-10, gap_extend=-0.5):
    """
    Calculate RMSDs for all given combinations
    """

    parser = PDBParser(QUIET=True)
    rows = []

    for comb in tqdm(combinations):

        pdb_a, chain_a, res_a, file_a = comb[0]
        pdb_b, chain_b, res_b, file_b = comb[1]

        ref_mod = parser.get_structure("X", str(file_a))
        mob_mod = parser.get_structure("Y", str(file_b))

        # Deal with NMR ensembles
        # Make all possible combinations of conformers in the files and compute
        # RMSDs on them
        ref_models = [ref_mod[idx][chain_a] for idx, _ in enumerate(ref_mod)]
        mob_models = [mob_mod[idx][chain_b] for idx, _ in enumerate(mob_mod)]
        model_combs = product(ref_models, mob_models)

        model_rmsds = []
        for model_comb in model_combs:
            reference, mobile = model_comb
            superimposition_rmsd = rmsd(reference, mobile, gap_open, gap_extend)
            model_rmsds.append(superimposition_rmsd)

        # If we have an NMR model, the final RMSD between the two files will be
        # the median of all RMSDs between conformers
        median_rmsd = np.median(model_rmsds)
        # rmsds.append(median_rmsd) #FIXME: can this be deleted?

        row = [pdb_a, chain_a, res_a,
               pdb_b, chain_b, res_b,
               median_rmsd, group]
        rows.append(row)

    return rows


def plot_histogram(df, out_path):
    """
    """
    sns.displot(data=df, x="RMSD", hue="Group", element="step")
    plt.xlabel("RMSD (Å)")
    sns.despine()
    dest_path = str(out_path / "rmsd_histogram.pdf")
    plt.savefig(dest_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_boxplot(df, out_path):
    """
    """
    ax = sns.boxplot(y="RMSD", x="Group", data=df,
                     notch=True, **{"bootstrap": 5000})
    ax = sns.swarmplot(y="RMSD", x="Group", data=df, color=".2", size=2)
    plt.ylabel("RMSD (Å)")
    dest_path = str(out_path / "rmsd_boxplot.pdf")
    plt.savefig(dest_path, dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":

    #########
    # Setup #
    #########

    args = digest_args()
    df = pd.read_csv(args.input_pairs)
    chains_path = args.structures
    num_workers = args.n_cpus
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

    ###################################
    # Compute all vs all global RMSDs #
    ###################################

    input_list = []
    print("Preparing input...")
    for phosphosite in psite_dict:
        between_combs = psite_dict[phosphosite]
        p_combs = phospho_combs[phosphosite]
        np_combs = nonphospho_combs[phosphosite]
        input_list.append(
            (phosphosite, between_combs, p_combs, np_combs, output_path))


    input_list = tuple(input_list)
    print("Ready to start RMSD calculations!")
    with Pool(num_workers) as pool:
        pool.starmap(compute_all_rmsds, tqdm(
            input_list, total=len(input_list)))
    pool.join()

    print("Finished!")
