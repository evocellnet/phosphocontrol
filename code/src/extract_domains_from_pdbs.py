"""
Script to extract domains from PDB-formatted structures
More concretely, this is used to extract the domains from the paired structures dataset.
"""

from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser

from utils.fasta import parse_fasta
from utils.structure import extract_sequence_from_chain, natural_idx_to_pdb_idx
from utils.alignment import align, map_seq_to_alignment, coverage
from utils.subs_mtx import BLOSUM62_LOOKUP
from utils.commands import call_selres

from tqdm import tqdm

from pprint import pprint
from warnings import warn


def make_structure_dict(subset_ptm_df, chains_path):
    """
    Create dictionary with protein structures and file paths
    {pdbcode_chain : ([Bio.PDB.Chain], file_path)}
    If the structure file contains multiple conformers (i.e. an NMR structure)
    the list will contain the different conformers
    """
    parser = PDBParser(QUIET=True)
    structure_dict = dict()
    
    for _, row in subset_ptm_df.iterrows():

        name_a = row["RENAMED_ENTITY_ONE"]
        _, chain_a = name_a.split("_")
        if name_a not in structure_dict.keys():
            a_file = chains_path / f"{name_a}.pdb"
            model_a = parser.get_structure("X", str(a_file))
            a_models = [model_a[idx][chain_a] for idx, _ in enumerate(model_a)]
            structure_dict[name_a] = (a_models, a_file)

        name_b = row["RENAMED_ENTITY_TWO"]
        _, chain_b = name_b.split("_")
        if name_b not in structure_dict.keys():
            b_file = chains_path / f"{name_b}.pdb"
            model_b = parser.get_structure("Y", str(b_file))
            b_models = [model_b[idx][chain_b] for idx, _ in enumerate(model_b)]
            structure_dict[name_b] = (b_models, b_file)

    return structure_dict


def get_structural_coverage(seq_aln, pdb_aln, seq_map, domain_start,
                            domain_end):
    """
    Calculate the structural coverage for a given region

    Arguments
    ---------
    seq_aln: aligned full sequence
    pdb_aln: aligned PDB sequence
    seq_map: dict, map of indexes in full sequence to alignment
    domain_start, domain_end: int, start and end positions of the domain

    Returns
    -------
    Float, structural coverage for the specified region
    """
    assert len(seq_aln) == len(pdb_aln)
    assert len(seq_aln) > 0

    try:
        mapped_domain_start = seq_map[domain_start]
        mapped_domain_end = seq_map[domain_end]
    except KeyError:
        print("KeyError here!")
        print("Domain start {domain_start}, end {domain_end}")
        print("Sequence alignment")
        pprint(seq_aln)
        print("PDB alignment")
        pprint(pdb_aln)
        print("Sequence map")
        pprint(seq_map)

    seq_domain = seq_aln[mapped_domain_start:mapped_domain_end]
    pdb_domain = pdb_aln[mapped_domain_start:mapped_domain_end]
    return coverage(seq_domain, pdb_domain)


def get_pdb_domain_idxs(domain_start, domain_end, seq_map, rev_pdb_map,
                        idx_to_pdb_idx):
    """
    Finds the PDB residue indexes corresponding to a specified range in the
    whole protein sequence

    Arguments
    ---------
    domain_start: int, start of the domain
    domain_end: int, end of the domain
    seq_map: mapping of sequence to alignment indexes
    rev_pdb_map: mapping of alignment to PDB sequence indexes (0-indexed)
    idx_to_pdb_idx: mapping of PDB sequence indexes (0-indexed)

    Returns
    -------
    pdb_domain_idxs: list of PDB residue indexes corresponding to the domain
    """
    domain_range = range(domain_start, domain_end)
    pdb_domain_idxs = []
    for idx in domain_range:
        try:
            # Get sequence position in the alignment
            aln_idx = seq_map[idx]
            # Position in the PDB alignment to PDB sequence
            pdb_aln_idx = rev_pdb_map[aln_idx]
            # Position in the PDB sequence to PDB residue index
            pdb_domain_idxs.append(idx_to_pdb_idx[pdb_aln_idx])
        except KeyError:
            pass
            # warn(f"Missing position {idx} in {v}", RuntimeWarning)

    return pdb_domain_idxs


if __name__ == "__main__":

    #########
    # Setup #
    #########

    base_path = Path("../../data/processed")

    # Structure pairs dataset
    ptm_df = pd.read_csv(base_path / "pdb_pairs" / "filtered_df.csv")
    # Mapping of phosphosites to Pfam domains
    psite_df = pd.read_csv(base_path / "phosphosite_to_pfamdomain" / "filtered_psite_per_domain.tsv")

    # Protein sequences
    fasta_f = base_path / "pdb_pairs" / "filtered_seqs.fasta"
    with open(fasta_f) as f:
        seqs = parse_fasta(f)
    # Structure input directory
    chains_path = base_path / "pdb_pairs" / "extracted_chains"/ "pdb"
    # Structure output directory
    output_path = base_path / "pdb_pairs" / "extracted_domains"
    output_path.mkdir(exist_ok=True)
    output_chains_path = base_path / "pdb_pairs" / "extracted_domains" / "pdb"
    output_chains_path.mkdir(exist_ok=True)

    ##############################
    # Retrieve domain structures #
    ##############################

    domain_coverage = []
    column_names = ["Uniprot ID", "Residue", "Phosphosite",
                    "PDB ID", "Chain ID", "Pfam ID", "Domain name", 
                    "Coverage", "Domain/sequence", "Domain/structure"]

    # Iterate over unique phosphosites
    for idx, row in tqdm(list(psite_df.iterrows())):
        psite = row["Phosphosite"]
        uniprot_id, res = psite.split("_")
        seq = seqs[uniprot_id]
        domain_id = row["Pfam ID"]
        domain_name = row["Domain name"]
        subset_ptm_df = ptm_df.loc[ptm_df["PHOSPHOSITE"] == psite]

        # Load structures pertaining to a given phosphosite into a dictionary
        # to avoid spending time on I/O
        structure_dict = make_structure_dict(subset_ptm_df,chains_path)

        for k, v in structure_dict.items():

            models, pdb_file = v[0], v[1]
            chain = models[0].get_id()

            # Pfam results are one-indexed, make them zero-indexed
            domain_start = int(row["Start"]) - 1
            domain_end = int(row["End"]) - 1
            domain_id = row["Pfam ID"]

            domain_len = len(list(range(domain_start, domain_end)))
            domain_seq_fraction = domain_len / len(seq)

            # Align Uniprot sequence to PDB sequence
            # Note: there should be no difference in sequence between different
            # conformers
            pdb_seq = extract_sequence_from_chain(models[0])
            seq_aln, pdb_aln = align(seq, pdb_seq, BLOSUM62_LOOKUP)
            seq_map = map_seq_to_alignment(seq, seq_aln)

            # Map PDB sequence positions to the alignment and vice versa
            pdb_seq_map = map_seq_to_alignment(pdb_seq, pdb_aln)
            rev_pdb_map = {v: k for k, v in pdb_seq_map.items()}

            # Map PDB sequence indexes to residue indexes
            idx_to_pdb_idx = natural_idx_to_pdb_idx(models[0])
            pdb_domain_idxs = get_pdb_domain_idxs(domain_start, domain_end,
                                                  seq_map, rev_pdb_map,
                                                  idx_to_pdb_idx)
            domain_structure_fraction = len(pdb_domain_idxs) / len(models[0])

            # Calculate coverage of the structure for the domain
            cov = get_structural_coverage(seq_aln, pdb_aln, seq_map,
                                          domain_start, domain_end)

            # Extract structure
            pdb_code, chain = pdb_file.stem.split("_")
            out_file = output_chains_path / f"{pdb_code}_{chain}_{domain_id}.pdb"
            # pdb_selres takes care of writing all conformers
            if len(pdb_domain_idxs) > 0:
                call_selres(str(pdb_file), str(out_file), pdb_domain_idxs)
            else:
                warn(f"No domain coverage for {pdb_code},{domain_id}", RuntimeWarning)

            # Add the data to the dataframe
            cov_row = [uniprot_id, res, psite,
                   pdb_code, chain, domain_id, domain_name,
                   cov, domain_seq_fraction, domain_structure_fraction]
            domain_coverage.append(cov_row)


    coverage_df = pd.DataFrame(domain_coverage, columns=column_names)
    df_path = str(output_path / "domain_coverage.csv")
    coverage_df.to_csv(df_path, index=False)

    print("All done here!")
