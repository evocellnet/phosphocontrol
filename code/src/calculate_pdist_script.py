"""
Script to perform pairwise calculations of residue-to-residue
distance matrix for
 - All phosphorylated versus non-phosphorylated structures
 - All phosphorylated versus non-phosphorylated structures
 - All non-phosphorylated versus non-phosphorylated structures
"""
import pandas as pd
import numpy as np
import warnings
import os
import sys
warnings.filterwarnings('ignore')

from scipy.spatial.distance import pdist

from numba import jit

from Bio import PDB

import psa.load as load
import psa.sequence as seq
import psa.elastic as elastic

# Functions
def calculate_pdist(name):
    xyz, label = load.coordinates(pps[name])
    xyz_p, label_p = pinfo[name]

    if label_p:
        label += label_p
        xyz = np.concatenate([xyz, xyz_p])

    idx = get_idx(label)

    dist = pdist(xyz)
    dist = np.concatenate([idx, dist])

    return dist


def get_phosphores(name, relabel=None, transform=None, model_id=0,
                   chain_ids=None, path='./', file_type='pdb'):
    
    # Load parser
    if file_type == 'pdb':
        parser = PDB.PDBParser()
    elif file_type == 'cif':
        parser = PDB.MMCIFParser()
    ppb = PDB.PPBuilder()
    
    # Download and load pdb/cif file
    load.download_file(name, path, file_type)
    filename = path + name + '.' + file_type
    structure = parser.get_structure(name, filename)

    # Load chains, maybe relabel
    chains = load.load_chains(structure[model_id], relabel, transform)
    
    # Select chain by id
    for ch in tuple(chains):
        if ch.id not in chain_ids:
            chains.remove(ch)
    
    coordinates = []
    labels = []

    # Loop over chains, peptides and residues
    for ch in chains:
        for res in ch:
                if res.get_resname() in ['SEP', 'TPO', 'PTR', 'HIP', 'NEP']:
                    coordinates.append(res['CA'].coord)
                    labels.append([res.full_id, res.resname])
    
    if len(coordinates) == 0:
        return np.array([]), []
    else:
        return np.array(coordinates, ndmin = 2), labels


def get_idx(labels):
    return np.array([label[0][3][1] for label in labels])

# argument list:

#01. psite ("uniprot"_"psite")
psite_arg = sys.argv[1]


# Load data
data = pd.read_csv("data/metadata/filtered_df.csv")


# Get list of proteins and PDBs
protP_accession = {}
protN_accession = {}

psites = list(set(data["PHOSPHOSITE"]))

for psite in psites:
    dataA = data[data["PHOSPHOSITE"] == psite]
    
    protein = psite.split("_")[0]
    
    protP_accession[psite] = list(set(dataA["AUTH_FULL_TWO"]))
    protN_accession[psite] = list(set(dataA["AUTH_FULL_ONE"]))
    
# Get list of pairwise comparisons

PP_list = {}
NN_list = {}
PN_list = {}

for group in data.groupby(["PHOSPHOSITE"]):
    df_group = group[1]
    
    psite = list(df_group["PHOSPHOSITE"])[0]
    P = list(df_group["AUTH_FULL_TWO"])
    N = list(df_group["AUTH_FULL_ONE"])
    
    PP_list[psite] = list(set(df_group["AUTH_FULL_TWO"]))
    NN_list[psite] = list(set(df_group["AUTH_FULL_ONE"]))
    PN_list[psite] = [i for i in zip(P, N)]

    
# Get PDB files

cif_list = ["6msb", "7ahz", "7ai1", "7ai0", "3ocb", "3qkm", "6l9t", 
            "8i9x", "8i9z", "8ia0", "8bp2", "7z9t", "7z6i", "7q5i", 
            "3bcr", "6yve", "8atl", "6weu", "6wew", "6wet", "6wfj", 
            "6wev", "7p0k", "7z0n", "7z3l", "7z3k", "7und", "7unc", 
            "7q6h", "8idt", "8idy", "8flc", "8inf", "8flb", "8fla", 
            "7oa0", "6yum", "6yul", "6z83", "6z84", "6tlu"]

pps = {}
pinfo = {}

psite = psite_arg
protein = psite.split("_")[0]
#index = psite.split("_")[1]
path = "data/PDB/" + protein + "/"

for group in [protP_accession, protN_accession]:
    for accession in group[psite]:
        name = accession.split("_")[0]
        chain = accession.split("_")[1]

        if name in cif_list:
            ft = "cif"
        else:
            ft = "pdb"


        try:
            structure = load.single_structure(name = name, 
                                              path = path,
                                              file_type = ft,
                                              chain_ids = "*",
                                              relabel = {chain: "*"})
            pps[accession] = structure

            xyz, labels = get_phosphores(name = name, 
                                         path = path,
                                         file_type = ft,
                                         chain_ids = "*",
                                         relabel = {chain: "*"})
            pinfo[accession] = [xyz, labels]

        except:
            group[psite].remove(accession)


# Save results
outpath = "output_pdist/"

if not os.path.exists(outpath):
    os.mkdir(outpath)

for pp in pps:
    dist = calculate_pdist(pp)
    with open(outpath + pp + "_pdist.npy", "wb") as f:
        np.save(f, dist)
