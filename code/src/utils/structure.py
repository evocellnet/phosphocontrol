"""
Module containing functions for dealing with Bio.PDB objects, create contact
matrices, and map contact matrices to multiple sequence alignments
"""

from math import sqrt
from dataclasses import dataclass
from pathlib import Path
import numpy as np

from utils.mappings import THREE_TO_ONE, LEGAL_AAS_THREE, ATOM_MASS
from utils.alignment import coverage

from warnings import warn


@dataclass(eq=True, frozen=True)
class StructureTuple:
    """
    Class for storing information about structures
    """

    pdb_code: str
    chain_id: str
    pdb_path: Path
    phospho_idxs: tuple  # Index of phosphosites in the alignment

#############################
# Bio.PDB object processing #
#############################

def remove_water(chain):
    """
    Given a protein chain, remove waters forming part of the sequence.

    Arguments
    ---------
    chain:  Bio.PDB.Chain.Chain object

    Returns
    -------
    filtered_chain: list of Bio.PDB.Residue.Residue, with no "water residues"
    """
    filtered_chain = []
    for res in list(chain):
        if res.resname == 'HOH':
            pass
        else:
            filtered_chain.append(res)
    return filtered_chain


def remove_non_legal_residues(chain, legal_aas=LEGAL_AAS_THREE):
    """
    Given a protein chain, remove non-legal amino acids forming part of the
    sequence.
    This is intended to remove non-amino acids from the chain.

    Arguments
    ---------
    chain:  Bio.PDB.Chain.Chain object

    Returns
    -------
    filtered_chain: list of Bio.PDB.Residue.Residue
    """
    filtered_chain = []
    for res in list(chain):
        if res.resname in legal_aas:
            filtered_chain.append(res)
    return filtered_chain


def extract_sequence_from_chain(chain, code=THREE_TO_ONE):
    """
    Translate a chain of amino acids in three letter code to one letter code.
    Note that the returned object contains no structural information: it is
    only the sequence of the protein.

    Called before "three_letter_to_one_letter"

    Arguments
    ---------
    chain: Bio.PDB.Chain.Chain object or list of Bio.PDB.Residue.Residue
    code:  dict, conversion table

    Returns
    -------
    translated_chain: string
    """
    translated_chain = []
    for res in list(chain):
        try:
            translated_chain.append(code[res.resname])
        except KeyError:
            # Ignore residues not in the code dictionary
            continue
    return ''.join(translated_chain)


def natural_idx_to_pdb_idx(chain):
    """
    Create a mapping from natural (zero-based) indexes to PDB chain indexes
    """
    return {idx: res.get_id()[1] for idx, res in enumerate(chain)}


def map_seq_pos_to_pdb_pos(pos, seq_map, rev_pdb_map, idx_to_pdb_idx, fallback):
    """
    Maps a given sequence index to the PDB residue number (according to the
    authors of the structure), after the sequence and the sequence derived from
    the PDB file have been aligned

    Arguments
    ---------
    pos: sequence position
    seq_map: mapping, sequence index -> alignment index
    rev_pdb_map: mapping, PDB sequence alignment index -> PDB sequence index
    idx_to_pdb_idx: mapping, PDB sequence index to PDB residue index
    fallback: hacky solution if residues are missing

    Returns
    -------
    pdb_pos: PDB residue number
    """
    # Map from position in the sequence to position in the alignment
    aln_pos = seq_map[pos]
    # Map to PDB sequence index
    try:
        pdbseq_pos = rev_pdb_map[aln_pos]
        # Map to PDB structure index
        pdb_pos = idx_to_pdb_idx[pdbseq_pos]
    except KeyError:
        try:
            if fallback == "forward":
                pdbseq_pos = rev_pdb_map[aln_pos + 1]
            elif fallback == "backward":
                pdbseq_pos = rev_pdb_map[aln_pos - 1]
            pdb_pos = idx_to_pdb_idx[pdbseq_pos]
        except KeyError:
            raise KeyError(f"Missing residue {pos} and neighbour")

    return pdb_pos


def structural_coverage(seq_aln, pdb_aln, seq_map, domain_start,
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
        raise ValueError(f"Missing position {domain_start} or {domain_end}")

    seq_domain = seq_aln[mapped_domain_start:mapped_domain_end]
    pdb_domain = pdb_aln[mapped_domain_start:mapped_domain_end]
    return coverage(seq_domain, pdb_domain)

#######################
# Abnormal structures #
#######################

def get_alpha_carbon_fraction(chain):
    """
    Get the fraction of residues in the chain that have an alpha carbon
    """
    ca_counter = 0
    skipped = 0
    for res in chain:
        resname = res.resname
        if resname in THREE_TO_ONE.keys():
            has_alpha = has_ca(res)
            if has_alpha:
                ca_counter += 1
        else:
            skipped += 1
    return ca_counter / (len(chain) - skipped)    

def has_ca(res):
    try:
        result = res["CA"]
        return True
    except KeyError:
        return False


def find_abnormal_residues(chain):
    abnormal_residues = 0
    skipped = 0
    for res in chain:
        resname = res.resname
        if resname in THREE_TO_ONE.keys():
            no_atoms = len(res)
            if no_atoms == 1:
                abnormal_residues += 1
        else:
            skipped += 1
    return abnormal_residues / (len(chain) - skipped)



###########################
# Structural descriptors  #
###########################

def radius_of_gyration(chain, atom_masses=ATOM_MASS):
    """
    Calculates the radius of gyration of a protein structure.

    Arguments
    ---------
    chain: Bio.PDB.Chain object

    Returns
    -------
    rg: float, radius of gyration
    """
    
    coords = []
    masses = []
    
    for res in chain:
        for atom in res:
            try:
                coords.append(list(atom.coord))
            except AttributeError:
                warn(f'Missing coordinates in residue {res._id[1]}')
            try:
                masses.append(atom_masses[atom.element])
            except KeyError:
                warn(f"Unknown element {atom.element} in residue {res._id[1]}")
    
    coord_mass_product = [(m*x, m*y, m*z) for (x, y, z), m in zip(coords, masses)]
    total_mass = sum(masses)
    rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coords, coord_mass_product))
    mm = sum((sum(i) / total_mass)**2 for i in zip(*coord_mass_product))
    rg = sqrt(rr / total_mass-mm)
    
    return rg


###########################
# Contact matrix building #
###########################

def calc_residue_distance(residue_one, residue_two):
    """
    Given two residues, calculate the Euclidean distance between them.
    The distance is measured between the beta carbons (or, in the case of
    glycine, with respect to the alpha carbon).

    Arguments
    ---------
    residue_one: Bio.PDB.Residue.Residue object
    residue_two: Bio.PDB.Residue.Residue object

    Returns
    -------
    euclid_dist:  float, Euclidean distance between the two residues
    """
    is_one_glycine = (residue_one.__dict__['resname'] == 'GLY')
    is_two_glycine = (residue_two.__dict__['resname'] == 'GLY')

    try:
        if is_one_glycine and is_two_glycine:
            diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
        elif is_one_glycine and not is_two_glycine:
            diff_vector = residue_one["CA"].coord - residue_two["CB"].coord
        elif not is_one_glycine and is_two_glycine:
            diff_vector = residue_one["CB"].coord - residue_two["CA"].coord
        else:
            diff_vector = residue_one["CB"].coord - residue_two["CB"].coord
    except KeyError: # Missing atom
        euclid_dist = np.nan
        return euclid_dist

    euclid_dist = np.sqrt(np.sum(diff_vector * diff_vector))
    return euclid_dist


def calc_intra_dist_matrix(chain):
    """
    Given one protein chain from a PDB file, calculate a matrix containing all
    pairwise Euclidean distances between residues within the chain.

    Arguments
    ---------
    chain: Bio.PDB.Chain.Chain object

    Returns
    -------
    dist_mtx: array-like, matrix containing pairwise Euclidean distances between
              residues, of dimensions (chain, chain)
    """

    dist_mtx = np.zeros((len(chain), len(chain)))
    for i, residue_a in enumerate(chain):
        for j, residue_b in enumerate(chain):
            if i != j: # Skip measuring distances to self
                dist_mtx[i, j] = calc_residue_distance(residue_a, residue_b)
    return dist_mtx


def calc_inter_dist_matrix(chain_one, chain_two):
    """
    Given two proteins chains from a PDB file, calculate a matrix containing
    all pairwise Euclidean distances between residues from the two different
    chains.

    Arguments
    ---------
    chain_one, chain_two:  Bio.PDB.Chain.Chain object, or
                           list of Bio.PDB.Residue.Residue

    Returns
    -------
    dist_mtx:   array-like, matrix containing pairwise Euclidean distances
                between residues, of dimensions (chain_one, chain_two)
    """
    dist_mtx = np.zeros((len(chain_one), len(chain_two)))
    for i, residue_one in enumerate(chain_one):
        for j, residue_two in enumerate(chain_two):
            dist_mtx[i, j] = calc_residue_distance(residue_one, residue_two)
    return dist_mtx


def get_contact_matrix(dist_mtx, threshold):
    """
    Discretize distance matrix according to a given threshold.
    This is done with the purpose of creating a contact map.

    Arguments
    ---------
    dist_mtx:   array-like, matrix containing pairwise Euclidean distances
                between residues, of dimensions (chain_one, chain_two)
    Returns
    -------
    contact_mtx: array-like, Boolean matrix; entries with a True value
                 represent a contact
    """
    return dist_mtx <= threshold

##################
# Matrix merging #
##################

def array_union(*args):
    """
    Given an arbitrary number of binary 2-dimensional arrays, find the
    union array.
    """

    # Check that all arrays are of the same dimension
    for array in args[1:]:
        assert args[0].shape == array.shape
    # Fill the array
    union_ar = np.zeros_like(args[0][0])
    for array in args:
        for idx, val in np.ndenumerate(array):
            if val == 1:
                idx = (idx[1], idx[2])
                union_ar[idx] = 1
    return union_ar


###################################
# Map contact matrix to alignment #
###################################

def map_contact_mtx_to_alignment(contact_mtx, aln_seq_a, aln_seq_b, mmap_a,
                                 mmap_b):
    """
    """

    # Get from the contact matrix which residues are in contact
    # for both sequences
    mapped_coords = map_contacts_to_alignment(
        contact_mtx, mmap_a, mmap_b)
    # Build a contact matrix adjusted to the dimensions of the aligned proteins
    mapped_contact_mtx = reconstruct_contact_mtx(
        mapped_coords, aln_seq_a, aln_seq_b)

    return mapped_contact_mtx


def map_contacts_to_alignment(contact_mtx, mmap_a, mmap_b):
    """
    Given a contact matrix, find what positions in the aligned sequences
    correspond to the contacts.

    Arguments
    ---------
    contact_mtx: array-like, contact matrix derived from the PDB file, of
                 dimensions (orig_seq_one, orig_seq_two)
    seq_map:     OrderedDict, in the format
                 {(index in original sequence, index in aligned sequence)}

    Returns
    -------
    mapped_coords: list of tuples, where each tuple (i,j) indicates a pair
                   of positions that make contact
    """
    contact_coords = get_contacts_from_matrix(contact_mtx)
    mapped_coords = []
    for contact in contact_coords:
        print(f"Contact: {contact}")
        try:
            pos_one = mmap_a[contact[0]]
            pos_two = mmap_b[contact[1]]
            if pos_one == "-" or pos_two == "-":
                print(f"""Contact {contact} contains deleted positions
                    (mapped as {(pos_one, pos_two)})""")
            else:
                mapped_coords.append((pos_one, pos_two))
        except KeyError:
            raise ValueError(f"Contact {contact} does not map to alignments")
    return mapped_coords


def map_seq_to_alignment(orig_seq, aln_seq):
    """
    Function to identify which positions in a sequence correspond to positions
    in the same sequence once it is aligned.

    Arguments
    ---------
    orig_seq: str, unaligned sequence
    aln_seq:  str, aligned sequence

    Returns
    -------
    seq_map: dict, in the format
             {(index in original sequence, index in aligned sequence)}

    """

    # Check that the two sequences are indeed the same,
    # and that no modifications have taken place during alignment
    # (e.g. MAFFT might delete positions when adding sequences to an alignment
    # if the --keeplength option is not passed)
    assert orig_seq == aln_seq.replace('-', '')

    seq_map = {}
    # Search each residue in the original sequence in the aligned sequence
    for i, res_one in enumerate(orig_seq):
        for j, res_two in enumerate(aln_seq):
            # They have to be the same residue and not have been mapped before
            if res_one == res_two and j not in seq_map.values():
                seq_map[i] = j
                break
    keys = seq_map.keys()
    vals = seq_map.values()

    # Check that the list of keys is sequential
    assert list(keys) == list(range(min(keys), max(keys) + 1))
    # Check that there are no duplicates in the values
    assert sorted(list(set(vals))) == list(vals)
    return seq_map


def get_contacts_from_matrix(contact_mtx):
    """
    Retrieve from the contact matrix the indexes at which there are contacts
    for the two sequences.

    Arguments
    ---------
    contact_mtx: array-like, contact matrix

    Returns
    -------
    contact_coords: list, contains tuples of indexes (i,j) indicating the
                    positions in the contact matrix where there are contacts
    """
    contact_coords = []
    for idx, val in np.ndenumerate(contact_mtx):
        if val == 1:
            contact_coords.append(idx)

    return contact_coords


def reconstruct_contact_mtx(mapped_coords, aln_seq_one, aln_seq_two):
    """
    Given indexes for the contacts corresponding to positions in the multiple
    sequence alignments, build a contact matrix of appropiate dimensions.

    Arguments
    ---------
    mapped_coords: list of tuples, where each tuple (i,j) indicates a pair
                   of positions that make contact
    aln_seq_one:    str, protein sequence orig_seq_one once it's been aligned
    aln_seq_two:    str, protein sequence orig_seq_two once it's been aligned

    Returns
    -------
    mapped_contact_mtx: array-like, contact matrix of dimensions
                        (aln_seq_one, aln_seq_two)
    """
    mapped_contact_mtx = np.zeros((len(aln_seq_one), len(aln_seq_two)))
    for contact in mapped_coords:
        mapped_contact_mtx[contact] = 1
    return mapped_contact_mtx


#################
# MAFFT mapping #
#################

def parse_mafft_map(mafft_map_path):
    """
    Read .map output of files, which specifies a mapping between positions
    in the original sequence and in the alignments. This includes positions
    (insertions with respect to the alignment) that have been deleted,
    which are indicated by the gap symbol.
    """
    with open(mafft_map_path) as f:
        data = f.readlines()
    data = [row.strip().split(",") for row in data[2:]]
    mmap = {}
    for row in data:
        if row[2] == " -":
            mmap[int(row[1])] = "-"
        else:
            mmap[int(row[1])] = int(row[2])
    return mmap
