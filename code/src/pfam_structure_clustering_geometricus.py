"""
Script to study where in the conformational landscape phosphorylated structures
are located.
"""

from pathlib import Path
import argparse
from random import seed
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

from utils import structure_embedding
from structure_embedding import compute_invariants, create_embedding
from warnings import warn

# Minimum total number of structures for a domain
global EXAMPLES_CUTOFF
EXAMPLES_CUTOFF = 10
# Minimum number of phosphorylated structures for a domain
global PHOSPHO_COUNT_CUTOFF
PHOSPHO_COUNT_CUTOFF = 4


def parse_args():
    """
    Parse script arguments
    """
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out_path", type=Path,
                    help="Output path")
    ap.add_argument("--seed", type=int, default=42,
                    help="RNG seed")
    args = ap.parse_args()
    return args


################
# Collect data #
################

def prepare_data(pfam_df, structures_path):
    """
    """
    # dict -> set of StructureTuple(pdb_code, chain_id, pdb_path,
    #                               tuple(phospho_pos)))
    pfam_to_structures = defaultdict(set)
    pfam_ids = list(pfam_df["Domain ID"].unique())
    for pfam_id in tqdm(pfam_ids):
        domain_path = structures_path / pfam_id
        subset_df = pfam_df.loc[pfam_df["Domain ID"]==pfam_id]
        for idx, row in subset_df.iterrows():
                pdb_id = row["PDB ID"]
                chain_id = row["Chain ID"]
                phospho_pos = str(row["Phosphosites"])
                pdb_path = str(domain_path / f"{pdb_id}_{chain_id}.pdb")

                if Path(pdb_path).is_file():
                    structure = structure_embedding.StructureTuple(pdb_id,chain_id,
                                                               pdb_path,phospho_pos)
                    pfam_to_structures[pfam_id].add(structure)

    return pfam_to_structures

if __name__ == "__main__":
    
    args = parse_args()
    out_path = args.out_path
    out_path.mkdir(exist_ok=True)
    seed(args.seed)

    base_path = Path("../../data/processed/pfam_structures/")
    structures_path = base_path / "extracted_domains_filtered"
    annot_df_f = base_path / "filtered_pfam_structures.tsv"
    annot_df = pd.read_csv(annot_df_f,delimiter='\t',index_col=0)
    subdirs = [x for x in structures_path.iterdir() if x.is_dir()]

    print("Preparing inputs...")
    pfam_to_structures = prepare_data(annot_df, structures_path)

    i = 0
    skipped_domains = []
    for subdir in subdirs:

        domain_name = subdir.stem
        print(f"Now on domain {domain_name}")

        domain_structures = list(pfam_to_structures[domain_name])

        domain_structures_dict = {}
        domain_structures = pfam_to_structures[domain_name]
        phospho_count = 0
        for domain_structure in domain_structures:
            pdb_code = domain_structure.pdb_code
            chain_id = domain_structure.chain_id
            phospho_idxs = domain_structure.phospho_idxs
            if phospho_idxs != "nan":
                phospho_count += 1
            domain_structures_dict[f"{pdb_code}_{chain_id}"] = domain_structure
        
        print(len(domain_structures), phospho_count)
        if len(domain_structures) >= EXAMPLES_CUTOFF and phospho_count >= PHOSPHO_COUNT_CUTOFF:

            domain_out_path = out_path / domain_name
            domain_out_path.mkdir(exist_ok=True)

            ###############################################
            # Cluster according to Geometricus embeddings #
            ###############################################

            print("Clustering according to Geometricus embedding...")

            invariants_kmer, invariants_radius = compute_invariants(domain_structures)
            embedding = create_embedding(invariants_kmer, invariants_radius)
            # embedding = pd.read_csv(domain_out_path / "geometricus_embedding.csv", index_col=0)
            print("Performing UMAP...")
            umap_reduced_embedding = structure_embedding.reduce_dimensionality_umap(embedding)
            umap_clusterer = structure_embedding.clustering(umap_reduced_embedding)
            umap_exemplars_dict = structure_embedding.get_cluster_exemplars(umap_clusterer)
            print("Performing PCA...")
            pca_reduced_embedding, pca = structure_embedding.reduce_dimensionality_pca(embedding)
            pca_clusterer = structure_embedding.clustering(pca_reduced_embedding)
            pca_exemplars_dict = structure_embedding.get_cluster_exemplars(pca_clusterer)

            #################
            # Create output #
            #################

            #structure_embedding.create_embedding_df(domain_structures, embedding, domain_out_path)
            structure_embedding.write_reduced_embedding(domain_structures, umap_reduced_embedding, domain_out_path, "umap")
            structure_embedding.write_reduced_embedding(domain_structures, pca_reduced_embedding, domain_out_path, "pca")

            structure_embedding.create_cluster_df(domain_structures, umap_clusterer, umap_exemplars_dict,
                                                  pca_clusterer, pca_exemplars_dict, domain_out_path)

            state_to_color = structure_embedding.assign_phospho_colors(domain_structures)

            structure_embedding.plot_clusters(umap_reduced_embedding, umap_clusterer, domain_name,
                                              domain_out_path, "UMAP")
            structure_embedding.plot_clusters(pca_reduced_embedding, pca_clusterer, domain_name,
                                              domain_out_path, "PCA",pca) 
            
            structure_embedding.plot_phospho_states(domain_structures, umap_reduced_embedding,
                                                    state_to_color, domain_out_path, domain_name, "UMAP")
            structure_embedding.plot_phospho_states(domain_structures, pca_reduced_embedding,
                                                    state_to_color, domain_out_path, domain_name, "PCA",
                                                    pca)

            i += 1
        
        else:
             warn(f"Domain {domain_name} has less than {EXAMPLES_CUTOFF} structures or {PHOSPHO_COUNT_CUTOFF} phosphorylated structures; skipping")
             skipped_domains.append(domain_name)
             continue

    with open(str(out_path / "skipped_domains.txt"), "w") as f:
        for domain_name in skipped_domains:
            to_write = "".join([domain_name, "\n"])
            f.write(to_write)

    print("Et voila!")
    print(f'Embedding carried out for {i} out of {len(subdirs)}')
