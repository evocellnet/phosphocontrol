# coding=utf-8
"""
Container for sets of sequences
"""

from utils.sequence import ProteinSequence
from utils.mappings import AA_TO_INT
from collections import defaultdict
from itertools import combinations
import numpy as np
from tqdm import tqdm

__author__ = "Miguel Correa"


class SequenceSet():
    """

    Attributes
    ----------
    seqs

    """

    def __init__(self):
        self.seqs = {}

    def __iter__(self):
        """Iterate over ProteinSequence objects"""
        for k, v in self.seqs.items():
            yield k, v

    def __len__(self):
        """Return number of sequences"""
        return len(self.seqs)

    def __bool__(self):
        return len(self) > 0

    def __setitem__(self, key, item):
        self.seqs[key] = item

    def __getitem__(self, key):
        return self.seqs[key]

    def __delitem__(self, key):
        del self.seqs[key]

    def keys(self):
        return self.seqs.keys()

    def values(self):
        return self.seqs.values()

    def items(self):
        return self.seqs.items()

    def pop(self, *args):
        return self.seqs.pop(*args)

    def __contains__(self, item):
        return item in self.seqs

    def __eq__(self, nonself):
        """
        """
        if isinstance(nonself, SequenceSet) and self.seqs == nonself.seqs:
            return True
        else:
            return False

    def add(self, sequence):
        """Adds a ProteinSequence object"""
        if isinstance(sequence, ProteinSequence):
            self.seqs[sequence.id] = sequence

    def read(self, lines):
        """Read FASTA file"""
        res = {}
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                label = line[1:]
                res[label] = []
            else:
                res[label].append(line)
        for k, v in res.items():
            seq = ProteinSequence(k, ''.join(v))
            self.add(seq)
        return self

    def is_alignment(self):
        """Verify whether all sequences are of the same length"""
        if len(self) > 1:
            lens = [len(v) for v in self.seqs.values()]
            if len(set(lens)) == 1:
                return True
            else:
                return False
        else:
            return False

    def aln_seq_size(self):
        if self.is_alignment:
            first = list(self.values())[0]
            return len(first)
        else:
            raise ValueError("Not an alignment")

    def del_gappy_cols(self, gap_cutoff):
        gap_num = AA_TO_INT['-']
        del_bool = []
        if self.is_alignment():
            num_mtx_t = self.num_matrix().T
            del_bool = []
            for idx, col in enumerate(num_mtx_t):
                if np.sum(col == gap_num) / len(col) > gap_cutoff:
                    del_bool.append(True)
                else:
                    del_bool.append(False)

            if np.all(del_bool):
                raise ValueError("All columns exceed the gap threshold")
            elif True in del_bool:
                self.del_cols(del_bool)
                return del_bool
            else:  # No gappy columns
                return []
        else:
            raise ValueError("Not an alignment")

    def del_constant_cols(self):
        if self.is_alignment():
            num_mtx = self.num_matrix()
            const_bool = np.bitwise_or.reduce(num_mtx) == num_mtx[0]
            if np.all(const_bool):
                raise ValueError("All columns are constant")
            elif True in const_bool:
                self.del_cols(const_bool)
                return const_bool
            else:  # No constant columns
                return self, []
        else:
            raise ValueError("Not an alignment")

    def del_cols(self, del_bool):
        """
        Delete a specified set of columns
        """
        for k, v in self.items():
            seq = "".join([aa for idx, aa in enumerate(v)
                           if del_bool[idx] == False])
            # TODO: this is probably slow
            seq = ProteinSequence(k, seq)
            self.add(seq)
        return self

    def num_matrix(self):
        seqs = list(self.seqs.values())
        no_seqs = len(seqs)
        len_seq = len(seqs[0])
        num_mtx = np.zeros((no_seqs, len_seq), dtype=np.int8)

        for idx, v in enumerate(seqs):
            num_mtx[idx] = v.seq_num()

        return num_mtx

    def bin_matrix(self):
        seqs = list(self.seqs.values())
        no_seqs = len(seqs)
        len_seq = len(seqs[0])
        bin_mtx = np.zeros((no_seqs, len_seq * 20), dtype=np.int8)

        for idx, v in enumerate(seqs):
            bin_mtx[idx] = v.seq_bin()

        return bin_mtx

    def num_dict(self):
        num_dict = {}
        for k, v in self.seqs.items():
            num_dict[k] = v.seq_num()
        return num_dict

    def bin_dict(self):
        bin_dict = {}
        for k, v in self.seqs.items():
            bin_dict[k] = v.seq_bin()
        return bin_dict

    def get_sample_weights(self, threshold=0.8):
        """
        Function to calculate weight of each sequence in the MSA.

        Arguments
        ---------
        threshold: float, identity above which sequences are clustered

        Returns
        -------
        sample_weights: dict, sample weights for the sequences in the alignment
                        {sequence id: sample weight}
        """
        if not self.is_alignment:
            raise ValueError("Not an alignment")

        if threshold < 0 or threshold > 1:
            raise ValueError(f"Identity threshold {threshold} out of bounds; should be within 0 and 1")

        combs = list(combinations(self.seqs.keys(), 2))
        idx_map = {k: idx for idx, k in enumerate(self.seqs.keys())}
        idents = defaultdict(lambda: defaultdict(int))

        for seq_i, seq_j in tqdm(combs):
            seq_id = self.seqs[seq_i].identity(self.seqs[seq_j])
            idents[seq_i][seq_j] = seq_id
            idents[seq_j][seq_i] = seq_id

        sample_weights = {}
        for k, idx in tqdm(idx_map.items()):
            seq_idents = np.array(list(idents[k].values()))
            try:
                w = 1 / \
                    np.where(seq_idents >= threshold)[0].size
            except ZeroDivisionError:
                w = 1
            sample_weights[k] = w

        return sample_weights

    def meff(self, threshold=0.8):
        """
        Calculate number of effective sequences
        """
        sw = self.get_sample_weights(threshold=threshold)
        return int(sum(sw))

    def column_aa_fractions(self, ignore_gaps=False):
        """
        Compute column-wise amino acid fraction.

        Output
        ------
        column_fractions: list of dicts containing amino acid fractions per
                          column, e.g,
                          [{0:0.25,1:0,...},..]

        """
        if not self.is_alignment:
            raise ValueError("Not an alignment")

        num_matrix = self.num_matrix()
        aa_fractions = []
        for row in num_matrix.T:
            column_fractions = {}
            unique, counts = np.unique(row, return_counts=True)
            counter = dict(zip(unique, counts))

            if ignore_gaps:
                try:
                    no_gaps = counter[20]
                except KeyError:
                    no_gaps = 0

            for k in counter.keys():
                if ignore_gaps and k != 20:
                    column_fractions[k] = counter[k] / \
                        (len(self) - no_gaps)
                elif not ignore_gaps:
                    column_fractions[k] = counter[k] / len(self)
            aa_fractions.append(column_fractions)
        return aa_fractions

    def column_entropy(self, ignore_gaps=True):
        """
        Calculate column-wise conservation in the alignment, understood as 
        Shannon information entropy.

        Arguments
        ---------
        ignore_gaps: whether to take into account gaps in the calculation

        Output
        ------
        column_entropy: array, Shannon entropy per column
        """
        column_entropy = []
        aa_fractions = self.column_aa_fractions(ignore_gaps=ignore_gaps)

        for col in aa_fractions:
            vals = list(col.values())
            column_entropy.append(-np.sum(vals * np.log(vals)))
        return column_entropy

