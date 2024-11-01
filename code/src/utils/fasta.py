"""
File containing utils for FASTA file I/O
"""

from warnings import warn

def parse_fasta(lines):
    """
    Return dict of {label:seq} from FASTA file

    Arguments
    ---------
    lines: list of lines, open file object

    Returns
    -------
    res: dict {label:seq}
    """
    res = {}  # result
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
        res[k] = ''.join(v)
    return res


def write_fasta(fasta_dict, out_path, flag='w'):
    """
    Write a FASTA file

    Arguments
    ---------
    fasta_dict: dict object, {header:sequence}
    out_path: str, path to output file
    flag: str, flag to open the file with

    Returns
    -------
    None

    """
    with open(out_path, flag) as f:
        for k, v in fasta_dict.items():
            f.write(''.join(['>', k, '\n']))
            f.write(''.join([v, '\n']))


def concatenate_fastas(file_list, out_path):
    """
    Read a list of FASTA files and concatenate them into one file.
    Makes sure that no headers are repeated.

    Arguments
    ---------
    file_list: list of FASTA files
    out_path: str, path to output file 

    Returns
    -------
    None
    """

    all_fasta = {}
    for file in file_list:
        with open(file) as f:
            fasta = parse_fasta(f)
            for header, seq in fasta.items():
                if header not in all_fasta.keys():
                    if len(seq) > 0:
                        all_fasta[header] = seq
                    else:
                        warn(f"Sequence {header} empty, skipping", RuntimeWarning)
    
    write_fasta(all_fasta, out_path)