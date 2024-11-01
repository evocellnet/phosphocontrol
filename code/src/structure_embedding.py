"""
Functions to create structural embeddings with Geometricus
"""
import numpy as np
import prody as prd

from geom.geometricus import MomentInvariants, SplitType, GeometricusEmbedding
from tqdm import tqdm
import umap


def compute_invariants(structures,kmer_split_size=16, radius_split_size=10):

    invariants_kmer = []
    invariants_radius = []

    for struct in tqdm(structures):

        key = (struct.pdb_code, struct.chain_id)
        pdb = prd.parsePDB(str(struct.pdb_path), chain=struct.chain_id)
        invariants_kmer.append(MomentInvariants.from_prody_atomgroup(
            key, pdb, split_type=SplitType.KMER, split_size=kmer_split_size))
        invariants_radius.append(MomentInvariants.from_prody_atomgroup(
            key, pdb, split_type=SplitType.RADIUS, split_size=radius_split_size))

    return invariants_kmer, invariants_radius


def create_embedding(invariants_kmer, invariants_radius, resolution=2.):
    kmer_embedder = GeometricusEmbedding.from_invariants(
        invariants_kmer, resolution=resolution)
    radius_embedder = GeometricusEmbedding.from_invariants(
        invariants_radius, resolution=resolution)
    embedding = np.hstack((kmer_embedder.embedding, radius_embedder.embedding))

    return embedding


def reduce_geometricus_dimensionality(embedding, metric="cosine",
                                      n_components=2, random_state=42):
    reducer = umap.UMAP(metric="cosine", min_dist=0.0,
                        n_components=n_components, random_state=random_state)
    geom_reduced = reducer.fit_transform(embedding)
    return geom_reduced