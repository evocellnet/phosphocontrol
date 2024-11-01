"""
Script to perform pairwise calculations of strain magnitude and
correlation dimension at the residue level for: 
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

from numba import jit
from scipy.ndimage import gaussian_filter1d

import psa.load as load
import psa.sequence as seq
import psa.elastic as elastic

# Functions
def calculate_CD(p1, p2):
    com_res, dict_1, dict_2 = seq.pairwise_alignment(pps[p1], pps[p2])
        
    xyz_1, label_1 = load.coordinates(pps[p1], com_res, dict_1)
    xyz_2, label_2 = load.coordinates(pps[p2], com_res, dict_2)
        
    weights = elastic.compute_weights_fast([xyz_1, xyz_2],
                                           parameters = [12]) #16

    try:
        F = elastic.deformation_gradient_fast(weights, xyz_1, xyz_2)
        

        _, gam_n = elastic.lagrange_strain(F)
        stretches, _ = elastic.principal_stretches_from_g(gam_n)

        lambda3 = np.absolute(stretches[:, 2])**.5
        idx = get_idx(label_1)

        order = np.argsort(lambda3)[::-1]
        Ns = range(10, len(lambda3))
        cd = np.zeros(len(Ns))

        for i, N in enumerate(Ns):
            cd[i] = correlation_dimension(N, xyz_1[order], xyz_2[order])
        
        return cd, lambda3, idx
        
    except:
        return np.nan

@jit(nopython = True)
def C(N, r, xyz_sorted):
    heaviside_sum = 0
    
    for i in range(N):
        for j in range(i+1, N):
            dist = 0
            for k in range(3):
                dist += (xyz_sorted[i][k] - xyz_sorted[j][k])**2
            dist = dist ** .5
            heaviside_sum += int(r > dist)
            
    return heaviside_sum / (N * (N - 1))

@jit(nopython = True)
def diam(xyz_sorted):
    diam = 0
    
    N = len(xyz_sorted)
    for i in range(N):
        for j in range(i+1, N):
            dist = 0
            for k in range(3):
                dist += (xyz_sorted[i][k] - xyz_sorted[j][k])**2
            
            if dist > diam:
                diam = dist
    return diam ** .5

def smoothed_gradient(N, xyz_sorted):
    r_range = np.linspace(8, diam(xyz_sorted), 60)
    Cs = np.array([C(N, r, xyz_sorted) for r in r_range])
    
    diff = np.zeros(len(r_range) - 1)
    r_vals = np.zeros(len(r_range) - 1)    
    for i in range(len(Cs) - 1):
        rp = r_range[i+1]
        rm = r_range[i]
        
        r_vals[i] = (rp + rm) / 2
        diff[i] = (np.log(Cs[i+1]) - np.log(Cs[i])) / (rp - rm) * r_vals[i]
        
    r_int = np.linspace(r_vals[0], r_vals[-1], 1000)
    
    diff_int = np.interp(r_int, r_vals, diff)
    diff_smt = gaussian_filter1d(diff_int, 15)
    
    return diff_smt, r_int

def correlation_dimension(N, xyz1_sorted, xyz2_sorted):
    grad1, r = smoothed_gradient(N, xyz1_sorted)
    grad2, r = smoothed_gradient(N, xyz2_sorted)
    
    return (np.max(grad1) + np.max(grad2)) / 2

def get_idx(labels):
    return np.array([label[0][3][1] for label in labels])

def compare_diff(ab_list, downsample = False, size = 20):
    n = len(ab_list)
    
    if downsample and n > int(size**2):
        lambda3 = np.zeros(int(size**2), dtype = object)
        
        random = np.random.choice(n, int(size**2), replace = False)
        ab_list = [ab_list[i] for i in random]
        
    else:
        lambda3 = np.zeros(n, dtype = object)
    
    
    c = 0
    for p1, p2 in ab_list:
        lambda3[c] = calculate_CD(p1, p2)
        c += 1
        
    return lambda3

def compare_same(aa_list, downsample = False, size = 20):
    n = len(aa_list)
    if n == 1:
        return np.nan
    
    if downsample and n > size:
        lambda3 = np.zeros(int(size*(size-1)), dtype = object)
        
        random = np.random.choice(n, size, replace = False)
        aa_list = [aa_list[i] for i in random]
        
    else:
        lambda3 = np.zeros(int(n*(n-1)), dtype = object)
        
    c = 0
    for p1 in aa_list:
        for p2 in aa_list:
            if p1 != p2:
                lambda3[c] = calculate_CD(p1, p2)
                c += 1
        
    return lambda3



# argument list:

#01. psite ("uniprot"_"psite")
psite_arg = sys.argv[1]

#02. downsample 0 = False, 1 = True
downsample_arg = sys.argv[2] == "1"

#03. downsampling size
size_arg = int(sys.argv[3])


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

        except:
            group[psite].remove(accession)


# Calculate strain

compPP = compare_same(PP_list[psite], downsample_arg, size_arg)
compNN = compare_same(NN_list[psite], downsample_arg, size_arg)
compPN = compare_diff(PN_list[psite], downsample_arg, size_arg)

# Save results
outpath = "output_cd_downsample_size_" + str(size_arg) + "/"
if not os.path.exists(outpath):
    os.mkdir(outpath)

with open(outpath + psite_arg + "_PP.npy", "wb") as f:
    np.save(f, compPP)
with open(outpath + psite_arg + "_NN.npy", "wb") as f:
    np.save(f, compNN)
with open(outpath + psite_arg + "_PN.npy", "wb") as f:
    np.save(f, compPN)

