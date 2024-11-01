# coding=utf-8
"""
Protein sequence container
"""

import numpy as np
from utils.mappings import CHR_MAPPING, AA_TO_INT
from Bio import pairwise2
from utils.subs_mtx import BLOSUM62_LOOKUP

__author__ = "Miguel Correa"


class ProteinSequence():
    """
    A class for storing and manipulating protein sequences. The object consists
    of an identifier and the sequence itself.

    The class is designed with contact prediction in mind. As contact
    prediction methods take into account the 20 standard proteinogenic amino
    acids plus the gap symbol, these are the characters supported in this
    class. Selenocysteine, pyrrolysine and ambiguous amino acids codes are not
    allowed.

    Attributes
    ----------
    id :  str, unique identifier
    seq : str, protein sequence
    """

    def __init__(self, ident, seq):
        """Initialize a ProteinSequence object"""
        self.id = ident
        self.seq = seq

    def __len__(self):
        """Sequence length"""
        return len(self.seq)

    def __bool__(self):
        return len(self) > 0

    def __iter__(self):
        for aa in self.seq:
            yield aa

    def __eq__(self, nonself):
        """
        Two ProteinSequence instances are equal iff they have the same id
        and the same sequence
        """
        if isinstance(nonself, ProteinSequence):
            is_equal = (self.id == nonself.id) and (self.seq == nonself.seq)
            return is_equal
        else:
            return False

    def __add__(self, nonself):
        """
        Concatenate ProteinSequence instances, creating a new instance
        """
        if isinstance(nonself, ProteinSequence):
            ident = '-'.join([self.id, nonself.id])
            seq = ''.join([self.seq, nonself.seq])
            return ProteinSequence(ident, seq)
        else:
            raise Exception('Other object is not a ProteinSequence instance')

    def __contains__(self, subseq):
        return subseq in self.seq

    def __getitem__(self, slice):
        return ProteinSequence(self.id, self.seq[slice])

    def __repr__(self):
        return self.seq

    def __str__(self):
        return self.seq

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, seq):
        # Check whether the sequences contain only allowed characters
        # (20 standard proteinogenic amino acids + gap symbol)
        if seq.translate(CHR_MAPPING) == seq:
            self._seq = seq
        else:
            raise ValueError('Unallowed amino acids in sequence')

    def replace(self, chr_a, chr_b, inplace=False):
        """Replace one string with another"""
        assert isinstance(chr_a, str)
        assert isinstance(chr_b, str)
        replaced = self.seq.replace(chr_a, chr_b)
        
        if inplace:
            self.seq = replaced
        else:
            return replaced

    def seq_ascii(self):
        """The protein sequence encoded in ASCII"""
        return bytearray(self.seq, "ascii")

    def seq_num(self):
        """The protein sequence as a series of ints."""
        return np.array([AA_TO_INT[aa] for aa in self.seq], dtype=np.int8)

    def seq_bin(self):
        """One-hot encoded protein sequence"""
        bin_seq = np.zeros(len(self.seq) * 20, dtype=np.int8)
        offset = 0
        for aa in self.seq:
            if aa == "-":  # Gaps are represented as 20 zeroes
                bin_seq[AA_TO_INT[aa] + offset] = 1
                offset += 20
        return bin_seq

    def hamming_distance(self, nonself):
        """
        Calculate Hamming distance (i.e. number of positions at which symbols
        are different between two strings of equal length) between two given
        ProteinSequences.
        """
        if isinstance(nonself, ProteinSequence):
            if len(self) == len(nonself):
                diff = np.bitwise_xor(self.seq_ascii(), nonself.seq_ascii())
                return np.nonzero(diff)[0].size
            else:
                raise Exception(
                    'Hamming distance or identity is only defined for strings of equal length')
        else:
            raise Exception('Other object is not a ProteinSequence instance')

    def identity(self, nonself):
        """
        Calculate identity between two given ProteinSequences.
        """
        hamming = self.hamming_distance(nonself)
        return 1 - (hamming / len(self))

    def similarity(self, nonself, sub_mtx=BLOSUM62_LOOKUP):
        """
        Calculate similarity between two given ProteinSequences.
        """
        if isinstance(nonself, ProteinSequence):
            match = 0
            # Check whether sequences are aligned
            if (len(self) == len(nonself)) and ("-" in self.seq) and ("-" in nonself.seq):
                for idx, aa in enumerate(self):
                    log_odds = score(aa, nonself.seq[idx], sub_mtx)
                    if log_odds > 0:
                        match += 1
                return match / len(self)
            else:
                seq1_aln, seq2_aln = self.align(nonself, sub_mtx)
                for idx, aa in enumerate(seq1_aln):
                    log_odds = score(aa, seq2_aln[idx], sub_mtx)
                    if log_odds > 0:
                        match += 1
                return match / len(seq1_aln)
        else:
            raise ValueError("Other object is not a ProteinSequence instance")

    def align(self, nonself, sub_mtx):
        """
        """
        seq1 = self.seq.replace("-", "")
        seq2 = nonself.seq.replace("-", "")
        aln = pairwise2.align.globaldx(seq1, seq2, sub_mtx)
        seq1_aln = aln[-1][0]
        seq2_aln = aln[-1][1]
        return seq1_aln, seq2_aln


def score(res1, res2, sub_mtx):
    """
    Return similarity score from BLOSUM matrix for two residues
    res1: string, amino acid
    res2: string, amino acid
    """
    return sub_mtx[(res1, res2)]
