"""
File containing aminoacid mappings.
Takes into account residues phosphorylated in eukaryotes (SEP,TPO,PTR),
as well as non-standard aminoacids (selenomethionine, pyrrolysine)
"""
from collections import defaultdict
from frozendict import frozendict

# Three-letter code (as one would find in e.g. a PDB file) to one-letter code
# Convert phosphorylated / modified residues to their non-phospho one-letter code
global THREE_TO_ONE
THREE_TO_ONE = frozendict({'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                           'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
                           'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
                           'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                           'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                           'SER': 'S', 'THR': 'T', 'TRP': 'W',
                           'TYR': 'Y', 'VAL': 'V', 'SEP': 'S',
                           'TPO': 'T', 'PTR': 'Y', 'MSE': 'M',
                           'PYL': 'K', 'SEC': 'C', '6BR': 'T',
                           'SVA': 'S', 'DTH': 'T', 'HIP': 'H',
                           'NEP': 'H', 'PHD': 'D'})

# Mapping of residues to an integer
global AA_TO_INT
AA_TO_INT = frozendict({'A': 0, 'R': 1, 'N': 2,
                        'D': 3, 'C': 4,
                        'E': 5, 'Q': 6, 'G': 7,
                        'H': 8, 'I': 9, 'L': 10,
                        'K': 11, 'M': 12, 'F': 13,
                        'P': 14, 'S': 15, 'T': 16,
                        'W': 17, 'Y': 18, 'V': 19,
                        '-': 20})


global INT_TO_AA
INT_TO_AA = frozendict({v:k for k,v in AA_TO_INT.items()})


global THREE_TO_INT
THREE_TO_INT = frozendict(
    {k: AA_TO_INT[THREE_TO_ONE[k]] for k in THREE_TO_ONE.keys()})

#######################
# Valid residue names #
#######################

global LEGAL_AAS_ONE
LEGAL_AAS_ONE = frozenset(['A', 'R', 'N',
                           'D', 'C', '-',
                           'E', 'Q', 'G',
                           'H', 'I', 'L',
                           'K', 'M', 'F',
                           'P', 'S', 'T',
                           'W', 'Y', 'V'])

global LEGAL_AAS_THREE
LEGAL_AAS_THREE = frozenset(THREE_TO_ONE.keys())

# For fast detection of non-legal residues
global CHR_MAPPING
CHR_MAPPING = defaultdict(str)
for aa in LEGAL_AAS_ONE:
    CHR_MAPPING[ord(aa)] = aa

global PHOSPHO_RESIDUES_THREE
PHOSPHO_RESIDUES_THREE = frozenset(["TPO", "SEP", "PTR","NEP","HIP"])

global PHOSPHO_RESIDUES_ONE
PHOSPHO_RESIDUES_ONE = frozenset(["T", "S", "Y", "H"])


###############
# Atom masses #
###############

global ATOM_MASS
ATOM_MASS = frozendict({'C': 12.0107,
                        'H': 1.0078,
                        'O': 15.9994,
                        'N': 14.0067,
                        'S': 32.065,
                        'P': 30.974,
                        'SE': 78.96})

########################
# Residue accesibility #
########################

# Maximum theoretical ASA in squared Angstroms values derived by :
# Tien et al 2013, Maximum Allowed Solvent Accessibilites of Residues in Proteins
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635#pone.0080635.s001

global RESIDUE_SURFACE
RESIDUE_SURFACES = frozendict({"A": 138,
                               "R": 285,
                               "D": 204,
                               "N": 204,
                               "C": 169,
                               "E": 233,
                               "Q": 234,
                               "G": 114,
                               "H": 231,
                               "I": 208,
                               "L": 211,
                               "K": 246,
                               "M": 227,
                               "F": 251,
                               "P": 166,
                               "S": 161,
                               "T": 182,
                               "W": 295,
                               "Y": 274,
                               "V": 184})
