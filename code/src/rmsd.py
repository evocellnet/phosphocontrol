"""
Sequence-based superimposition of two protein structures

Note: direct use of this script is not setup to deal with NMR structures properly - 
it only uses the first model of the ensemble.

all_global_rmsds (using the functions here) does take NMR structures into account
"""

import argparse
from pathlib import Path

from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio import pairwise2

from utils.mappings import THREE_TO_ONE
from utils.subs_mtx import BLOSUM62_LOOKUP

from warnings import warn

def digest_args():
    """
    Parse script arguments
    """

    ap = argparse.ArgumentParser(description=__doc__)

    ap.add_argument("--ref", type=Path, help="Reference structure")
    ap.add_argument("--mob", type=Path, help="Mobile structure")
    ap.add_argument("--r_chain", type=str, help="Reference structure chain")
    ap.add_argument("--m_chain", type=str, help="Mobile structure chain")
    ap.add_argument("-o", "--output", type=Path, help="Output path")
    ap.add_argument("--gap_open", type=int, help="Gap open penalty",
                    default=-10)
    ap.add_argument("--gap_extend", type=int, help="Gap extend penalty",
                    default=-0.5)

    args = ap.parse_args()
    return args


def get_indexed_sequence(chain, aa_map=THREE_TO_ONE):
    """

    Arguments
    ---------
    chain:  Bio.PDB.Chain object
    aa_map: dict, mapping between three-letter aminoacid code to one-letter

    Returns
    -------
    idxed_seq: list of tuples, [(residue idex, one letter residue)]

    """
    idxed_seq = [(res.id[1], aa_map[res.__dict__['resname']]) for res in chain if res.__dict__['resname'] in aa_map.keys()]
    return idxed_seq


def create_idx_chain(chain):
    """
    Function to workaround the fact that biopython does not allow you access
    HETATMs (which PTMs often are) by index.

    Arguments
    ---------
    chain:  Bio.PDB.Chain object

    Returns
    -------
    idx_chain: dict, {PDB index: Bio.PDB.Residue.Residue objects}
    """
    idx_chain = {res.get_id()[1]: res for res in chain}
    return idx_chain


def align_sequences(seq_a, seq_b, gap_open, gap_extend,
                    sub_mtx=BLOSUM62_LOOKUP):
    """
    Perform global pairwise alignment between two sequences
    """

    assert (len(seq_a) > 0) and (len(seq_b) > 0)
    aln = pairwise2.align.globalds(seq_a, seq_b, sub_mtx, gap_open, gap_extend,
                                   penalize_end_gaps=(False, False),
                                   one_alignment_only=True)
    seq_a_aln = aln[-1][0]
    seq_b_aln = aln[-1][1]
    return seq_a_aln, seq_b_aln


def map_sequences(idx_seq_a, idx_seq_b, seq_a_aln, seq_b_aln):
    """
    Create mapping between the two protein chains
    """
    mapping = {}
    aa_i_A, aa_i_B = 0, 0

    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(seq_a_aln, seq_b_aln)):

        if aa_aln_A == "-":
            if aa_aln_B != "-":
                aa_i_B += 1
        elif aa_aln_B == "-":
            if aa_aln_A != "-":
                aa_i_A += 1
        else:
            assert idx_seq_a[aa_i_A][1] == aa_aln_A
            assert idx_seq_b[aa_i_B][1] == aa_aln_B
            mapping[idx_seq_a[aa_i_A][0]] = idx_seq_b[aa_i_B][0]
            aa_i_A += 1
            aa_i_B += 1

    return mapping

def rmsd(reference, mobile, gap_open, gap_extend):

    idx_ref_seq = get_indexed_sequence(reference)
    idx_mob_seq = get_indexed_sequence(mobile)

    ref_seq = "".join([x[1] for x in idx_ref_seq])
    mob_seq = "".join([x[1] for x in idx_mob_seq])

    ref_aln, mob_aln = align_sequences(ref_seq, mob_seq, gap_open,
                                       gap_extend)

    mapping = map_sequences(idx_ref_seq, idx_mob_seq, ref_aln, mob_aln)

    # Work around biopython
    ref_idx_chain = create_idx_chain(reference)
    mob_idx_chain = create_idx_chain(mobile)

    # Create lists of matching C-alphas
    ref_ca_list, mob_ca_list = [], []
    for ref_res, mob_res in mapping.items():

        try:
            ref_ca = ref_idx_chain[ref_res]["CA"]
            mob_ca = mob_idx_chain[mob_res]["CA"]

            ref_ca_list.append(ref_ca)
            mob_ca_list.append(mob_ca)

        except KeyError:
            # Alpha carbons might be missing from the structure
            #warn(f"Missing alpha carbon at {ref_res,mob_res}", RuntimeWarning)
            pass

    # Superimpose matching residues
    if not ref_ca_list or not mob_ca_list:
        # No matching atoms
        warn(f"No matching atoms; returning NaN")
        return float("nan")
    else:
        si = Superimposer()
        si.set_atoms(ref_ca_list, mob_ca_list)
        si.apply(mobile.get_atoms())
        return si.rms


if __name__ == "__main__":

    args = digest_args()
    # Parse structures & take only the necessary chain
    p = PDBParser(QUIET=True)
    struct_reference = p.get_structure("A", str(args.ref))

    if len(struct_reference) > 0:
        warn("Reference structure is an ensemble; using only the first model")

    try:
        reference = struct_reference[0][args.r_chain]
    except KeyError:
        raise Exception(f"Chain {args.r_chain} not found in reference \
                          structure")

    struct_mobile = p.get_structure("B", str(args.mob))

    if len(struct_mobile) > 0:
        warn("Mobile structure is an ensemble; using only the first model")

    try:
        mobile = struct_mobile[0][args.m_chain]
    except KeyError:
        raise Exception(f"Chain {args.m_chain} not found in mobile structure")

    # Align sequences and get mapping between residues
    rms = rmsd(reference, mobile, args.gap_open, args.gap_extend)

    print(f"RMSD between structures: {rms}")
