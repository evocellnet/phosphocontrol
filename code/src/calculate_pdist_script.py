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

########################
#                      #
# Function Definitions #
#                      #
########################


def calculate_pdist(name):
    """
    Compute the pairwise distances between alpha-carbon atoms in a protein structure.
    If phosphorylated residues are present, include their coordinates as well.

    Parameters:
        name (str): Unique identifier for the protein chain (e.g., PDBID_PSITE).

    Returns:
        np.ndarray: Vector containing the residue indices followed by pairwise 
                    distances between C-alpha atoms.
    """
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
    """
    Extract coordinates of phosphorylated residues (e.g., SEP, TPO, 
    PTR) from a PDB or CIF file.

    Parameters:
        name (str): PDB ID
        relabel (dict): Optional chain relabeling
        transform (dict): Optional transformation
        model_id (int): Model number (usually 0)
        chain_ids (list/str): Specific chain(s) to extract
        path (str): Path to the structure file
        file_type (str): 'pdb' or 'cif'

    Returns:
        tuple: (coordinates (np.ndarray), labels (list))
    """
    
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

    # Search for phosphorylated or modified residues
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
    """
    Extract residue indices from label metadata.

    Parameters:
        labels (list): Metadata list for residues

    Returns:
        np.ndarray: Residue indices
    """
    return np.array([label[0][3][1] for label in labels])

#########################
#                       #
# Main Script Execution #
#                       #
#########################

# Command line arguments
psite_arg = sys.argv[1]                      # e.g., "P12345_67"

# Load phosphosite metadata
data = pd.read_csv("data/metadata/filtered_df.csv")


# Organize P/N structures per phosphosite
protP_accession = {}
protN_accession = {}

psites = list(set(data["PHOSPHOSITE"]))

for psite in psites:
    dataA = data[data["PHOSPHOSITE"] == psite]
    
    protein = psite.split("_")[0]
    
    protP_accession[psite] = list(set(dataA["AUTH_FULL_TWO"]))
    protN_accession[psite] = list(set(dataA["AUTH_FULL_ONE"]))
    
# Define list of structures available in CIF format
cif_list = ["6msb", "7ahz", "7ai1", "7ai0", "3ocb", "3qkm", "6l9t", 
            "8i9x", "8i9z", "8ia0", "8bp2", "7z9t", "7z6i", "7q5i", 
            "3bcr", "6yve", "8atl", "6weu", "6wew", "6wet", "6wfj", 
            "6wev", "7p0k", "7z0n", "7z3l", "7z3k", "7und", "7unc", 
            "7q6h", "8idt", "8idy", "8flc", "8inf", "8flb", "8fla", 
            "7oa0", "6yum", "6yul", "6z83", "6z84", "6tlu"]

# Load structures and phosphorylated residue information

pps = {}
pinfo = {}

psite = psite_arg
protein = psite.split("_")[0]
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
            # Load structure using custom loader
            structure = load.single_structure(name = name, 
                                              path = path,
                                              file_type = ft,
                                              chain_ids = "*",
                                              relabel = {chain: "*"})
            pps[accession] = structure
            
            # Extract phosphorylated residues            
            xyz, labels = get_phosphores(name = name, 
                                         path = path,
                                         file_type = ft,
                                         chain_ids = "*",
                                         relabel = {chain: "*"})
            pinfo[accession] = [xyz, labels]

        except:
            # If structure fails, remove from list
            group[psite].remove(accession)


# Save results
outpath = "output_pdist/"

if not os.path.exists(outpath):
    os.mkdir(outpath)

for pp in pps:
    dist = calculate_pdist(pp)
    with open(outpath + pp + "_pdist.npy", "wb") as f:
        np.save(f, dist)

