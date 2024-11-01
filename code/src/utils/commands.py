"""
Contains functions to call CLI programs
"""

import os
import subprocess
from Bio.PDB.DSSP import DSSP
import pathlib
from pathlib import Path

from warnings import warn

######################
# STRUCTURE MAPPINGS #
######################

global SS_DICT
SS_DICT = {"H": "Alpha helix", "B": "Isolated beta-bridge residue",
           "E": "Strand", "G": "3-10 helix", "I": "Pi helix",
           "T": "Turn", "S": "Bend", "-": "Loop/irregular"}


######################
# Sequence alignment #
######################

def call_muscle(in_file, out_file):

    call = ["muscle", "-in", str(in_file), '-out', str(out_file)]
    try:
        subprocess.check_call(call, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as error:
        raise ValueError(f'MUSCLE failed: {error}')
    except OSError:
        raise OSError("MUSCLE executable not found in PATH")
    

def call_clustalo(in_file, out_file, hmm_in, distmat_out=None, threads=2,
                  outfmt='fasta'):
    if distmat_out:
        call = ["clustalo", "--dealign", "--force", "--in", in_file,
                "--hmm-in", hmm_in, "--out", out_file, "--outfmt", outfmt,
                "--full", "--distmat-out", distmat_out,
                "--threads", str(threads), "-v"]
    else:
        call = ["clustalo", "--dealign", "--force", "--in", in_file,
                "--hmm-in", hmm_in, "--out", out_file, "--outfmt", outfmt,
                "--full", "--threads", str(threads), "-v"]

    try:
        subprocess.check_call(call, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as error:
        raise ValueError(f'clustalo failed: {error}')
    except OSError:
        raise OSError("clustalo executable not found in PATH")


def call_clustalo_nohmm(in_file, out_file, distmat_out=None, threads=2,
                        outfmt='fasta'):
    if distmat_out:
        call = ["clustalo", "--dealign", "--force", "--in", in_file,
                "--out", out_file, "--outfmt", outfmt,
                "--full", "--distmat-out", distmat_out,
                "--threads", str(threads), "-v"]
    else:
        call = ["clustalo", "--dealign", "--force", "--in", in_file,
                "--out", out_file, "--outfmt", outfmt,
                "--full", "--threads", str(threads), "-v"]

    try:
        subprocess.check_call(call, stdout=subprocess.DEVNULL)
        # subprocess.check_call(call)
    except subprocess.CalledProcessError as error:
        raise ValueError(f'clustalo failed: {error}')
    except OSError:
        raise OSError("clustalo executable not found in PATH")


###################################
# Structural alignment/similarity #
###################################


def call_tmalign(pdb_a, pdb_b, out_file):
    """
    Function to call TMalign; captures stdout and redirects it to a file.
    """
    call = ["TMalign", pdb_a, pdb_b]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"TMalign failed: {error}")
    except OSError:
        raise OSError("TMalign binary not found in PATH")


def call_caretta(input_pdb_folder, out_path, threads=8, full=False, 
                 matrix=False):
    """
    Writes alignment features
    """

    call = ["caretta-cli", input_pdb_folder, "-o",
            out_path, "-t", str(threads), "--features"]
    if full:
        call.append('--full')
    if matrix:
        call.append('--matrix')

    try:
        subprocess.check_call(call)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"caretta-cli failed: {error}")
    except OSError:
        raise OSError("caretta-cli binary not found in PATH")

def process_caretta_alignment(filename):
    """
    Remove the .pdb ending from caretta fasta headers
    """

    call = ['sed','-i','-e','s/.pdb//g' ,filename]

    try:
        subprocess.check_call(call)
    except subprocess.CalledProcessError as error:
        raise ValueError(f'sed failed: {error}')
    except OSError:
        raise OSError('sed or file not found')

######################
# Structure commands #
######################

def call_dssp(model, pdb_path, dssp_path, ss_dict=SS_DICT):
    """

    Arguments
    ---------
    model: Bio.PDB.Model.Model
    pdb_path: path to PDB file
    dssp_path: path to DSSP binary
    ss_dict: secondary structure mapping from DSSP output

    Return
    ------
    dssp_data

    """
    dssp = DSSP(model, pdb_path, dssp=dssp_path)
    dssp_data = {}
    keys = list(dssp.keys())
    for k in keys:
        data = dssp[k]
        idx, aa, ss, asa = data[0], data[1], data[2], data[3]
        ss = ss_dict[ss]
        dssp_data[f"{aa}{idx}"] = [ss, asa]

    return dssp_data

#############
# pdb-tools #
#############


def pdb_combo(in_dir, out_dir):
    """
    Call pdb_call_tidy_res for each PDB file in a directory.
    """

    pdbs = list(in_dir.glob("*.pdb"))
    # Double check that we do have the expected input
    assert len(pdbs) > 0

    for pdb in pdbs:

        pdb_id = pdb.stem
        out_file = out_dir / f"{pdb_id}.pdb"
        pdb_call_tidy_reres(str(pdb), str(out_file))

def pdb_call_tidy_reres(in_file, out_file):
    """
    Call a combo a pdbtools utilities to get around the bug with Caretta where
    unusual modified aminoacids are not 
    """
    #sort_call = ["pdb_sort", in_file]
    tidy_call = ["pdb_tidy", in_file]
    reres_call = ["pdb_reres", "-0"]

    with open(out_file, 'w') as outf:

        # process_sort = subprocess.Popen(sort_call, stdout=subprocess.PIPE,
        #                                 shell=False)
        process_tidy = subprocess.Popen(tidy_call,
                                        stdout=subprocess.PIPE, shell=False)
        process_reres = subprocess.Popen(reres_call, stdin=process_tidy.stdout,
                                         stdout=outf, shell=False)
        process_tidy.stdout.close()

    # pdb_sort 4CXU_A.pdb | pdb_tidy | pdb_reres -0 > test.pdb

def call_pdbtofasta(in_file, out_file):
    """
    Function to call pdb_tofasta (from pdb-tools)
    """

    call = ["pdb_tofasta", str(in_file)]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f,stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"pdb_tofasta failed: {error}")
    except OSError:
        raise OSError("pdb_tofasta binary not found in PATH")


def call_pdbfetch(pdb_code, out_file, biounit=False):
    """
    Function to call pdb_fetch (from pdb-tools)
    """

    assert len(pdb_code) == 4
    warn("Asserting PDB ID has 4 characters; will have 8 in the future",
         DeprecationWarning)
    if biounit:
        call = ["pdb_fetch", pdb_code]
    else:
        call = ["pdb_fetch", "-biounit", pdb_code]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"pdb_fetch failed: {error}")
    except OSError:
        raise OSError("pdb_fetch binary not found in PATH")


def call_selchain(in_file, chains, out_file):
    """
    Function to call pdb_selchain (from pdb-tools)

    Arguments
    ---------
    pdb: path to input file (.pdb, .cif)
    chains: list, chain identifiers
    out_file: path to output file

    Returns
    -------
    None; writes file to disk
    """

    chains = ",".join(chains)
    chains = "".join(["-", chains])
    call = ["pdb_selchain", chains, in_file]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f,stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"pdb_selchain failed: {error}")
    except OSError:
        raise OSError("pdb_selchain binary not found in PATH")


def call_reres(in_file, out_file, start):
    """
    Call pdb_reres to renumber the residues of a PDB file starting from a
    given number (can be negative)
    """

    assert isinstance(start, int)

    # print("call_reres here")
    # print(in_file)
    # print(out_file)

    in_file = str(pathlib.Path.expanduser(Path(in_file)))
    out_file = str(pathlib.Path.expanduser(Path(out_file)))

    call = ["pdb_reres", f"-{start}", in_file]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"pdb_reres failed: {error}")
    except OSError:
        raise OSError("pdb_reres binary not found in PATH")


def call_selres(in_file, out_file, res_list):
    """
    Function to call pdb_selres (from pdb-tools) to select residues by index,
    provided as res_list.

    This function does not offer the full flexibility of calling pdb_selres
    (e.g. does not allow you specify multiple ranges)

    pdb_selres deals properly with NMR structures,
    i.e., it uses all models of the ensemble, not just the first one
    """
    assert len(res_list) > 0
    res_list = ",".join([str(x) for x in res_list])
    res_list = "".join(["-", res_list])
    call = ["pdb_selres", res_list, in_file]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"pdb_selres failed: {error}")
    except OSError:
        raise OSError("pdb_selres binary not found in PATH")


def call_selres_range(in_file, out_file, res_start, res_end):
    """
    Function to call pdb_selres to select a range of residues
    """
    call = ["pdb_selres", f"-{res_start}:{res_end}", in_file]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"pdb_selres failed: {error}")
    except OSError:
        raise OSError("pdb_selres binary not found in PATH")


def create_ensemble(in_files, out_file):
    call_mkensemble(in_files, out_file)
    is_ensemble_ok = check_ensemble(out_file)
    if not is_ensemble_ok:
        raise ValueError(f"Error in PDB ensemble {out_file}")


def call_mkensemble(in_files, out_file):
    """
    Create ensemble of structures using pdb_mkensemble
    """
    assert len(in_files) > 1

    in_files = " ".join([in_files])
    call = ["pdb_mkensemble", in_files]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"pdb_mkensemble failed: {error}")
    except OSError:
        raise OSError("pdb_mkensemble binary not found in PATH")


def check_ensemble(in_file):
    """
    Check the validity of a file containing an ensemble using pdb_chkensemble
    """
    call = ["pdb_chkensemble", in_file]

    try:
        subprocess.check_call(
            call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as error:
        return False
    except OSError:
        raise OSError("pdb_chkensemble binary not found in PATH")

    return True


def call_oneinall(infile, chains, out_file):
    """
    In one go:
        - Extracts specified chains
        - Tidy up the PDB file

    Uses shell = True!

    pdb_delhetatm will delete phosphorylated residues!

    Arguments
    ---------
    pdb: path to input file (pdb or cif)
    chains: list, chain identifiers
    """
    chains = ",".join(chains)
    call = f"pdb_selchain -{chains} {infile} | pdb_tidy > {out_file}"

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f, shell=True)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"call_oneinall failed: {error}")
    except OSError:
        raise OSError("pdb_tools binaries not found in PATH")


#########
# hmmer #
#########


def call_hmmalign(in_fasta, hmm_file, out_file):
    """
    """

    call = ["hmmalign", "--amino", hmm_file, in_fasta]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"hmmalign failed: {error}")
    except OSError:
        raise OSError("hmmalign binary not found in PATH")


def call_hmmfetch(pfam_db, domain_name, out_file):
    """
    """
    call = ["hmmfetch", pfam_db, domain_name]

    try:
        with open(out_file, "w") as f:
            subprocess.check_call(call, stdout=f)
    except subprocess.CalledProcessError as error:
        raise ValueError(f"hmmfetch failed: {error}")
    except OSError:
        raise OSError("hmmfetch binary not found in PATH")


def call_hmmsearch(in_file, hmm_model, out_prefix, nobias=False, n_cpus=4,
                   seed=42):
    """
    Call hmmsearch

    Arguments
    ---------

    Returns
    -------
    None, writes files to disk
    """
    if nobias:
        call = ["hmmsearch", "--seed", str(seed), "--cpu",
                str(n_cpus), "--nobias", "-A", f"{out_prefix}.sto"
                "--tblout", f"{out_prefix}.tsv",
                hmm_model, in_file]
    else:
        call = ["hmmsearch", "--seed", str(seed), "--cpu",
                str(n_cpus), "-A", f"{out_prefix}.sto",
                "--tblout", f"{out_prefix}.tsv",
                hmm_model, in_file]
    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(call, stdout=devnull)
    except subprocess.CalledProcessError:
        raise ValueError('hmmsearch failed')
    except OSError:
        raise OSError('hmmsearch executable not found in PATH')


def call_esl_reformat(in_file, out_file, seq_format="fasta"):
    """
    Call esl-reformat to turn Stockholm files into FASTA files
    """

    call = ["esl-reformat", seq_format, in_file]
    try:
        with open(out_file, "w") as seqfile:
            subprocess.check_call(call, stdout=seqfile)
    except subprocess.CalledProcessError:
        raise ValueError('esl-reformat failed')
    except OSError:
        raise OSError('esl-reformat executable not found in PATH')

###########
# hhsuite #
###########


def call_hhblits(fasta_path, uniclust_path, evalue_cutoff, output_path):
    """
    Intended to generate features for Cfold
    """

    call = ["hhblits", "-i", str(fasta_path), "-d", str(uniclust_path),
            "-E", str(evalue_cutoff), "-cpu", "1", "-all", "-oa3m",
            str(output_path)]
    try:
        subprocess.check_call(call)
    except subprocess.CalledProcessError:
        raise ValueError('hhblits failed')
    except OSError:
        raise OSError('hhblits executable not found in PATH')

