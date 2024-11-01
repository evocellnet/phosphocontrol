"""
Functions to create structural embeddings with Geometricus
"""
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass

from geometricus import get_invariants_for_structures, Geometricus, SplitInfo, SplitType
import umap
import hdbscan
import pandas as pd
import numpy as np

from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns

from threadpoolctl import threadpool_limits


@dataclass(eq=True, frozen=True)
class StructureTuple:
    """
    Class for storing information about structures
    """
    pdb_code: str
    chain_id: str
    pdb_path: str
    # Index(es) of phosphosite(s)
    phospho_idxs: str


#####################
# Create embeddings #
#####################

def embed_structures(structures, kmer_split_size=16,
                     radius_split_size=8, resolution=6., n_threads=8):
    """
    Creates a Geometricus embedding for a set of PDB files.

    Arguments
    ---------
    pdb_files: list of StructureTuple objects
    kmer_split_size: int, kmer chunk size
    radius_split_size: int, radius around fragment
    resolution: float, resolution at which to discretize invariants
    n_threads: int, number of threads to use for embedding

    Returns
    -------
    count_matrix: np array, shapemer embedding
    """
    names = [f"{struct.pdb_code}.pdb" for struct in structures]
    pdb_files = [Path(struct.pdb_path).resolve() for struct in structures]

    kmer_matrix = get_kmer_invariants(
        pdb_files, names, kmer_split_size, resolution, n_threads)
    radius_matrix = get_radius_invariants(
         pdb_files, names, radius_split_size, resolution, n_threads)

    # Create matrix with all shapemers
    #count_matrix = radius_matrix
    count_matrix = np.hstack([kmer_matrix, radius_matrix])
    return count_matrix


def get_kmer_invariants(pdb_files, names, kmer_split_size, 
                        resolution, n_threads):
    """
    """
    invariants, _ = get_invariants_for_structures(pdb_files, n_threads=n_threads,
                                                  split_infos=[
                                                      SplitInfo(SplitType.KMER, kmer_split_size)],
                                                  moment_types=["O_3", "O_4", "O_5", "F"])

    shapemer_class = Geometricus.from_invariants(
        invariants, resolution=resolution)

    _, shapemers = shapemer_class.map_protein_to_shapemer_indices()
    kmer_matrix = shapemer_class.get_count_matrix(shapemer_keys=shapemers)

    return kmer_matrix


def get_radius_invariants(pdb_files, names, radius_split_size, resolution, n_threads):
    """
    """
    invariants, _ = get_invariants_for_structures(pdb_files, n_threads=n_threads,
                                                  split_infos=[
                                                      SplitInfo(SplitType.RADIUS, radius_split_size)],
                                                  moment_types=["O_3", "O_4", "O_5", "F"])

    shapemer_class = Geometricus.from_invariants(
        invariants, resolution=resolution)
    _, shapemers = shapemer_class.map_protein_to_shapemer_indices()
    radius_matrix = shapemer_class.get_count_matrix(shapemer_keys=shapemers)

    return radius_matrix

############################
# Dimensionality reduction #
############################

def reduce_dimensionality_umap(embedding, metric="cosine",
                                      n_components=2, random_state=42):
    reducer = umap.UMAP(metric=metric, min_dist=0.0,
                        n_components=n_components, random_state=random_state)
    reduced_embedding = reducer.fit_transform(embedding)
    return reduced_embedding

def reduce_dimensionality_pca(embedding, n_components=2, random_state=42):
    """
    """
    # Workaround for bug where PCA gets stuck 
    with threadpool_limits(limits=1):
        pca = PCA(n_components, random_state=random_state)
        reduced = pca.fit_transform(embedding)
        return reduced, pca



##############
# Clustering #
##############

def clustering(X,metric="euclidean"):
    """
    Perform clustering with HDBSCAN using Euclidean distance
    """

    if X.shape[0] > 10:
        clusterer = hdbscan.HDBSCAN(metric=metric, prediction_data=True,
                                    min_samples=10)
    else: # Few samples
        clusterer = hdbscan.HDBSCAN(metric=metric, prediction_data=True,
                                    min_samples=5)

    clusterer.fit(X)
    return clusterer


def get_cluster_exemplars(clusterer):
    """
    Get the indexes of samples that HDBSCAN consider as exemplars of each
    cluster.

    Returns
    -------
    exemplars: dict, {cluster: [sample_idxs]}

    """
    selected_clusters = clusterer.condensed_tree_._select_clusters()
    raw_condensed_tree = clusterer.condensed_tree_._raw_tree

    exemplars = []
    for cluster in selected_clusters:

        cluster_exemplars = np.array([], dtype=np.int64)
        for leaf in clusterer._prediction_data._recurse_leaf_dfs(cluster):
            leaf_max_lambda = raw_condensed_tree['lambda_val'][
                raw_condensed_tree['parent'] == leaf].max()
            points = raw_condensed_tree['child'][
                (raw_condensed_tree['parent'] == leaf) &
                (raw_condensed_tree['lambda_val'] == leaf_max_lambda)]
            cluster_exemplars = np.hstack([cluster_exemplars, points])
        exemplars.append(cluster_exemplars)

    exemplars_dict = {}
    for idx, cluster_exemplars in enumerate(exemplars):
        exemplars_dict[idx] = cluster_exemplars.tolist()

    return exemplars_dict

##########
# Output #
##########

def create_embedding_df(domain_structures, embedding, out_path):
    """
    """
    row_names = [f"{struct.pdb_code}_{struct.chain_id}" for struct in domain_structures]
    # Write Geometricus embedding
    embedding_df = pd.DataFrame(embedding, index=row_names)
    embedding_df.to_csv(out_path / "geometricus_embedding.csv")


def create_cluster_df(domain_structures, umap_clusterer, umap_exemplars_dict,
                      pca_clusterer, pca_exemplars_dict, out_path):
    """
    """
    # TODO: would be good to add uniprot IDs
    columns = ["PDB ID","Chain ID","Phosphosites","UMAP cluster","UMAP is exemplar",
               "PCA cluster", "PCA is exemplar", "File path"]
    cluster_df = []

    for idx, struct in enumerate(domain_structures):

        umap_cluster = umap_clusterer.labels_[idx]
        pca_cluster = pca_clusterer.labels_[idx]
        umap_is_exemplar = check_if_exemplar(idx, umap_exemplars_dict, umap_cluster)
        pca_is_exemplar = check_if_exemplar(idx, pca_exemplars_dict, pca_cluster)
        
        pdb_path = str(Path(struct.pdb_path))

        row = [struct.pdb_code, struct.chain_id,
               struct.phospho_idxs, str(umap_cluster), 
               str(umap_is_exemplar), str(pca_cluster),
               str(pca_is_exemplar), pdb_path]
        cluster_df.append(row)
    cluster_df = pd.DataFrame(cluster_df,columns=columns)
    cluster_df.to_csv(out_path / "cluster_df.csv",sep='\t',index=None)


def check_if_exemplar(idx, exemplars_dict, cluster):
    # Noise points (-1) cannot be exemplar!
    if (cluster != -1) and (idx in exemplars_dict[cluster]):
        is_exemplar = True
    else:
        is_exemplar = False

    return is_exemplar


def write_reduced_embedding(domain_structures, reduced_embedding, out_path, method):
    
    row_names = [f"{struct.pdb_code}_{struct.chain_id}" for struct in domain_structures]
    columns = ["Component 1","Component 2"]
    reduced_embedding_df = pd.DataFrame(index=row_names)
    reduced_embedding_df["component_1"] = reduced_embedding[:,0]
    reduced_embedding_df["component_2"] = reduced_embedding[:,1]
    reduced_embedding_df.to_csv(out_path / f"{method}_reduced_geometricus_embedding.csv")


def create_embedding_df(domain_structures, embedding, out_path):
    """
    """
    row_names = [f"{struct.pdb_code}_{struct.chain_id}" for struct in domain_structures]
    # Write Geometricus embedding
    embedding_df = pd.DataFrame(embedding, index=row_names)
    embedding_df.to_csv(out_path / "geometricus_embedding.csv")


#########
# Plots # 
#########

def plot_clustermap(embedding, phospho_colors, out_path):
    """
    Plot clustermap using average linkage hierarchical clustering
    """
    sns.clustermap(embedding, method="average",
                   metric="cosine", xticklabels=False,
                   yticklabels=False,
                   row_colors=phospho_colors)

    plot_path = str(out_path / "embedding_clustermap.pdf")
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()



def plot_clusters(reduced_embedding, clusterer, domain_name, out_path, method, 
                  reducer_object=None):
    """
    Create scatterplot of UMAP embedding with cluster labels
    """
    sns.set_context('talk')
    sns.set_style('white')
    colors = ListedColormap(
        list(sns.color_palette("husl", len(set(clusterer.labels_)))))
    plt.figure(figsize=(9, 9))

    scatter = plt.scatter(
        reduced_embedding[:, 0], reduced_embedding[:, 1], alpha=0.3, 
        c=clusterer.labels_,
        cmap=colors)
    plt.legend(handles=scatter.legend_elements()[0],
               labels=set(clusterer.labels_),
               bbox_to_anchor=(1.05, 1))

    if method == 'UMAP':
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")
    elif method == "PCA":
        if reducer_object is not None:
            plt.xlabel(f"PC 1 ({100*reducer_object.explained_variance_ratio_[0]:.1f}%)")
            plt.ylabel(f"PC 2 ({100*reducer_object.explained_variance_ratio_[1]:.1f}%)")
        else:
            plt.xlabel("PC 1")
            plt.ylabel("PC 2")
    else:
        raise ValueError(f"Invalid method {method}")

    plt.title(f"{domain_name}")
    sns.despine()
    plt.savefig(out_path / f"{method}_{domain_name}_clusters.pdf",
                dpi=150,bbox_inches='tight')
    plt.close()


def plot_phospho_states(domain_structures, reduced_embedding,
                        state_to_color, out_path, domain_name, method,
                        reducer_object=None):
    """
    """
    sns.set_context('poster')
    sns.set_style('white')
    plt.figure(figsize=(10, 10))

    phospho_state_to_structure = defaultdict(set)
    for idx, domain_structure in enumerate(domain_structures):
        phospho_idxs = domain_structure.phospho_idxs
        if phospho_idxs == "nan":
            phospho_idxs = ""
        phospho_state_to_structure[phospho_idxs].add(idx)

    for phospho_state, idxs in phospho_state_to_structure.items():
        if phospho_state == "":
            alpha = 0.4
        else:
            alpha = 1

        x = [item for idx, item in enumerate(reduced_embedding[:, 0]) if idx in idxs]
        y = [item for idx, item in enumerate(reduced_embedding[:, 1]) if idx in idxs]

        sns.set_context("talk")
        plt.scatter(x, y, alpha=alpha, s=200,
                    c=np.array([state_to_color[phospho_state]]),
                    label=str(phospho_state))
    
    plt.legend(loc="best", title="Phosphorylation state",
               bbox_to_anchor=(1.05, 1),frameon=False,fontsize=24,
               prop={'size':22})

    if method == 'UMAP':
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")
    elif method == "PCA":
        if reducer_object is not None:
            plt.xlabel(f"PC 1 ({100*reducer_object.explained_variance_ratio_[0]:.1f}%)",fontsize=30)
            plt.ylabel(f"PC 2 ({100*reducer_object.explained_variance_ratio_[1]:.1f}%)",fontsize=30)
        else:
            plt.xlabel("PC 1")
            plt.ylabel("PC 2")
    else:
        raise ValueError(f"Invalid method {method}")

    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)


    plt.title(f"{domain_name}",fontsize=32)
    sns.despine()
    plt.savefig(out_path / f"{method}_{domain_name}_phospho_state.pdf",
                dpi=120,bbox_inches='tight')
    plt.close()


def assign_phospho_colors(domain_structures, palette_type="colorblind"):
    """
    Returns
    -------
    phospho_colors:             list, color for each domain structure
    phospho_state_to_structure: dict,
                                {(phospho state):[idx_struct_a,
                                                  idx_struct_b,...]}
    state_to_color:             dict,
                                {(phospho state):rgb color tuple}
    """

    # Determine the different phosphorylation states
    # Store the indexes of the structures in each one
    phospho_states = set()
    phospho_state_to_structure = defaultdict(list)
    structure_to_phospho_state = {}
    for idx, struct in enumerate(domain_structures):
        phospho_state = struct.phospho_idxs
        if phospho_state == "nan":
            phospho_state = ""
        phospho_states.add(phospho_state)
        #phospho_state_to_structure[phospho_state].append(idx)
        #structure_to_phospho_state[idx] = phospho_state
    phospho_states = sorted(list(phospho_states))

    # Assign a color to each phospho state
    palette = sns.color_palette(palette_type, len(phospho_states))
    state_to_color = {}
    for idx, state in enumerate(phospho_states):
        state_to_color[state] = palette[idx]

    # phospho_colors = []
    # for idx, struct in enumerate(domain_structures):
    #     state = structure_to_phospho_state[idx]
    #     phospho_colors.append(state_to_color[state])

    return state_to_color
    #return phospho_colors, phospho_state_to_structure, state_to_color
    #return phospho_state_to_structure, state_to_color
