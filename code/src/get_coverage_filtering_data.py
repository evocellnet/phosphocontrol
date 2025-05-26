#!/usr/bin/env python
"""
This script processes structural data and full-length UniProt sequences to calculate two
metrics for each structure chain in the dataset:
    1. Coverage: fraction of the full sequence that is represented in the structure.
    2. Median similarity: average sequence similarity of a structure to other structures
       of the same protein.
The output is used to filter out structures for downstream analyses.

Usage: python calculate_overlap.py

Note: this script calls Clustal Omega to create multiple sequence alignments.
      It must be installed in the system for it to run.

Inputs
------
    - ../../data/raw/pdb_pairs/phosphorylated_residues_filtered_quality.csv:
      CSV containing chain-to-UniProt mappings and PDB entity identifiers.
    
    - ../../data/raw/pdb_pairs/uniprot/seqs.fasta:
      Multi-FASTA file containing full-length UniProt sequences.

    - ../../data/raw/pdb_pairs/extracted_chains/pdb/*.pdb:
      Directory of PDB files for individual chains.

Outputs
-------
    - ../../data/raw/pdb_pairs/sequences/<UNIPROT_ID>/*.fasta:
      FASTA files of structure-derived and full sequences, per protein.

    - ../../data/raw/pdb_pairs/sequences/<UNIPROT_ID>/alignment.fasta:
      Multiple sequence alignments (MSA) of structure chains per protein.

    - ../../data/raw/pdb_pairs/overlap_df.csv:
      Final dataset containing coverage and similarity values for each chain.

"""

from pathlib import Path
import pandas as pd
import numpy as np
import Bio.PDB
from Bio.PDB import PDBParser

from utils.structure import extract_sequence_from_chain
from utils.alignment import coverage
from utils.sequenceset import SequenceSet
from utils.fasta import parse_fasta, write_fasta, concatenate_fastas
from utils.commands import call_clustalo_nohmm

from tqdm import tqdm
# import pickle


def get_overlap_per_chain(unique_proteins, sequences_path):
    """
    Calculate overlap metrics for each polypeptide chain for each protein

    Arguments
    ---------
    unique_proteins: list, list of Uniprot IDs
    sequences_path:  path to directory containing subdirs for each Uniprot ID
                     The function assumes that all sequences necessary have
                     already been written to their corresponding subdirs
    Returns
    -------
    overlap_df: pandas DataFrame
                entity_id -> PDB ID + chain ID
                coverage -> structural coverage with respect to the full sequence
                median_overlap -> median structural overlap with respect to every
                                  other structure for this sequence
    """
    overlap_df = []
    columns = ["uniprot_id","pdb_id","entity_id","coverage","median_similarity","sequence_length"]
    for unique_protein in tqdm(unique_proteins):

        if unique_protein != "P30305":
        
            print(f"Now aligning {unique_protein}")
            protein_path = sequences_path / unique_protein
            # Pool all sequences for the protein into one FASTA file
            fasta_files = list(protein_path.glob("*.fasta"))
            nonaln_seqs = str(protein_path / "all_sequences.fasta")
            aln_seqs = str(protein_path / "alignment.fasta")

            if not Path(aln_seqs).is_file():
                concatenate_fastas(fasta_files, nonaln_seqs)
                # Align sequences
                call_clustalo_nohmm(nonaln_seqs, aln_seqs, threads=8)
            else:
                print('Already aligned!')
            alignment = SequenceSet()
            with open(aln_seqs) as f:
                alignment.read(f)

            ref_seq = alignment[unique_protein]
            del alignment[unique_protein]
            similarity_dict = {}
            
            similarity_matrix = np.zeros((len(alignment),len(alignment)))
            for i, (header_i, seq_i) in enumerate(alignment.items()):
                if header_i != unique_protein:
                    # Get coverage wrt the full sequence
                    cov = coverage(str(ref_seq), str(seq_i))
                    seq_length = len(str(seq_i).replace("-",""))
                    
                    for j, (header_j, seq_j) in enumerate(alignment.items()):
                        if i == j:
                            similarity_matrix[i,j] = 1
                            similarity_matrix[j,i] = 1
                        elif header_j != unique_protein:

                            try:
                                # Check whether we have computed this one previously
                                similarity = similarity_dict[header_j][header_i]
                            except KeyError:
                                similarity = seq_i.similarity(seq_j)
                                similarity_dict[header_i] = {header_j:similarity}
                                similarity_matrix[i,j] = similarity
                                similarity_matrix[j,i] = similarity
                    
                    median_similarity = np.median(similarity_matrix[i,:])
                    #print(f"Median similarity for {header_i}: {median_similarity}")

                    pdb_id = header_i.split("_")[0]
                    row = [unique_protein, pdb_id, header_i, cov, median_similarity, seq_length]
                    overlap_df.append(row)
    
    print("Finished alignments!")
    overlap_df = pd.DataFrame(overlap_df, columns=columns)
    return overlap_df


if __name__ == "__main__":

    #########
    # Setup #
    #########

    data_path = Path("../../data/raw/pdb_pairs")
    uniprot_path = data_path / "uniprot"
    df = pd.read_csv(data_path / "phosphorylated_residues_filtered_quality.csv")
    chains_path = data_path / "extracted_chains" / "pdb"

    ################################
    # Write full protein sequences #
    ################################
    
    print("Writing Uniprot sequences...")
    # Get the full sequences from Uniprot
    with open(str(uniprot_path / "seqs.fasta")) as f:
        all_seqs = parse_fasta(f)
    # Keep only the Uniprot ID
    all_seqs = {k.split("|")[1]:v for k,v in all_seqs.items()}
    unique_proteins = list(df["UNIPROT"].unique())

    # Create a directory where to write protein sequences
    sequences_path = data_path / "sequences"
    sequences_path.mkdir(exist_ok=True)
    # Write each Uniprot sequence
    for unique_protein in unique_proteins:
        protein_path = sequences_path / unique_protein
        protein_path.mkdir(exist_ok=True)
        seq_dict = {}
        seq_dict[unique_protein] = all_seqs[unique_protein]
        write_fasta(seq_dict, str(protein_path / "full_sequence.fasta"))

    ###################################
    # Write sequences from structures #
    ###################################
    
    # Extract sequences from chains and write them to their corresponding directory
    df["AUTH_FULL_ONE"] = df["ENTRY_ID_ONE"] + "_" + df["RENAMED_ASYM_ID_ONE"]
    df["AUTH_FULL_TWO"] = df["ENTRY_ID_TWO"] + "_" + df["RENAMED_ASYM_ID_TWO"]

    set_one = set(list(df["AUTH_FULL_ONE"]))
    set_two = set(list(df["AUTH_FULL_TWO"]))
    full_set = list(set_one.union(set_two))

    # # Associate each chain with its corresponding UniProt ID
    chain_to_uniprot_id = {}
    print("Associating chains to Uniprot IDs...")
    for chain_name in tqdm(full_set):
        subset_df = df.loc[(df["AUTH_FULL_ONE"]==chain_name) | (df["AUTH_FULL_TWO"]==chain_name)]
        uniprot_id = list(subset_df["UNIPROT"].unique())[0]
        chain_to_uniprot_id[chain_name] = uniprot_id

    # For debugging
    # with open('chain_to_uniprot_id.pickle', 'wb') as handle:
    #     pickle.dump(chain_to_uniprot_id, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Extract sequences from chains    
    print("Extracting sequences from structures...")
    parser = PDBParser(QUIET=True)
    for chain_name in tqdm(full_set):
        chain_id = chain_name.split("_")[1]
        try:
            structure = parser.get_structure("0000", str(chains_path / f"{chain_name}.pdb"))
            chain = structure[0][chain_id]
            seq = extract_sequence_from_chain(chain)
            chain_dict = {}
            chain_dict[chain_name] = seq
            protein_path = sequences_path / chain_to_uniprot_id[chain_name]
            print(f"Writing {chain_name} to {protein_path}")
            write_fasta(chain_dict, str(protein_path / f"{chain_name}.fasta"))
        except Bio.PDB.PDBExceptions.PDBConstructionException:
            print(f"Error processing file for {chain_name}")

    #######################
    # Calculate coverages #
    #######################

    print("Calculating coverages...")
    overlap_df = get_overlap_per_chain(unique_proteins, sequences_path)
    overlap_df.to_csv(data_path / "overlap_df.csv",index=False)

                         

        



