"""
Contains functions related to protein sequence alignment
"""

from Bio import pairwise2


def align(seq_a, seq_b, sub_mtx, gap_open=-10, gap_extend=-0.5):
    """
    Align two protein sequences
    """
    assert (len(seq_a) > 0) and (len(seq_b) > 0)
    aln = pairwise2.align.globalds(seq_a, seq_b, sub_mtx, gap_open, gap_extend,
                                   penalize_end_gaps=(False, False),
                                   one_alignment_only=True)
    seq_a_aln = aln[-1][0]
    seq_b_aln = aln[-1][1]
    return seq_a_aln, seq_b_aln


def coverage(ref_seq, seq):
    """
    Calculate the coverage provided by a given (possibly partial) sequence with
    respect to a reference sequence, once they have been aligned
    """
    assert len(ref_seq) == len(seq)
    assert len(ref_seq) > 0
    matches = 0
    gaps = 0
    for idx, res in enumerate(ref_seq):
        if res != "-" and seq[idx] != "-":
            matches += 1
        # Account for possible gaps in the reference sequence
        # (e.g. insertions in PDB files with respect to Uniprot sequence)
        elif res == "-" and seq[idx] != "-":
            gaps += 1
    cov = matches / (len(ref_seq)-gaps)

    return cov


def similarity(seq1, seq2, sub_mtx):
    """
    Calculate the similarity between two aligned protein sequences
    """
    assert len(seq1) == len(seq2)
    match = 0
    for idx, res in enumerate(seq1):
        log_odds = score(res, seq2[idx], sub_mtx)
        if log_odds > 0:
            match += 1
    return match / len(seq1)


def score(res1, res2, sub_mtx):
    return sub_mtx[(res1, res2)]


def map_seq_to_alignment(orig_seq, aln_seq):
    """
    Function to identify which positions in a sequence correspond to positions
    in the same sequence once it is aligned.

    TODO: This could be made faster (no need for nested for loop)

    Arguments
    ---------
    orig_seq: str, unaligned sequence
    aln_seq:  str, aligned sequence

    Returns
    -------
    seq_map: dict, in the format
             {index in original sequence : index in aligned sequence}
    """

    seq_map = {}
    for i, res_one in enumerate(orig_seq):
        for j, res_two in enumerate(aln_seq):
            if res_one == res_two and j not in seq_map.values():
                seq_map[i] = j
                break
        else:
            pass

    # Sanity checks
    vals = seq_map.values()
    # Check that the list of values increases monotonically
    assert is_increasing_monotonic(list(vals))
    # Check that there are no duplicates in the values
    assert sorted(list(set(vals))) == list(vals)

    return seq_map


def is_increasing_monotonic(array):
    """
    Returns whether the the elements in a list increase monotonically or not
    """
    return all(x < y for x, y in zip(array, array[1:]))
