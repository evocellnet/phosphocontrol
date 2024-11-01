"""
Script to fetch the all the structures of Pfam domains found in the dataset.
To do so, it uses the PDB structures associated with the domain in SIFTS.
"""

from pathlib import Path
import pandas as pd
from time import sleep
from tqdm import tqdm
from warnings import warn
from Bio.PDB import PDBParser
import prody as prd
from collections import defaultdict

from utils.fasta import parse_fasta, write_fasta
from utils.commands import call_selres_range, call_oneinall, call_clustalo
from utils.structure import extract_sequence_from_chain, natural_idx_to_pdb_idx, map_seq_pos_to_pdb_pos
from utils.alignment import align, map_seq_to_alignment

#from utils.mappings import PHOSPHO_RESIDUES_THREE

from predict_status.helper import find_phospho_idxs

def query_pdb(pdb_id, out_path):
    """
    Retrieve a PDB file from the PDB
    
    Prody does not like Path objects, use strings!
    Also, it doesn't want you to provide a name for the file,
    just the name of the directory
    """

    assert len(pdb_id) == 4
    pdb_success = prd.fetchPDB(pdb_id, compressed=False, folder=str(out_path))
    sleep(0.1)
    return pdb_success


if __name__ == "__main__":

    #########
    # Setup #
    #########

    # Read file with phosphosite to Pfam associations
    psite_df = pd.read_csv(
        "../../data/processed/phosphosite_to_pfamdomain/filtered_psite_per_domain.tsv")
    unique_pfam_domains = list(psite_df["Pfam ID"].unique())

    ptm_df = pd.read_csv("../../data/processed/pdb_pairs/filtered_df.csv")

    sifts_path = Path("../../data/raw/sifts")
    pdb_pfam_mapping = pd.read_csv(sifts_path / "pdb_pfam_mapping.csv",skiprows=1,low_memory=False)
    pdb_pfam_mapping = pdb_pfam_mapping.loc[pdb_pfam_mapping["PFAM_ACCESSION"].isin(unique_pfam_domains)]
    
    pdb_chain_pfam = pd.read_csv(sifts_path / "pdb_chain_pfam.csv",skiprows=1,low_memory=False)
    pdb_chain_pfam = pdb_chain_pfam.loc[pdb_chain_pfam["PFAM_ID"].isin(unique_pfam_domains)]

    # Directories for output
    out_path = Path("../../data/processed/pfam_structures")
    out_path.mkdir(exist_ok=True)
    full_path = out_path / "full_structures"
    full_path.mkdir(exist_ok=True)
    chains_path = out_path / "extracted_chains"
    chains_path.mkdir(exist_ok=True)
    domains_path = out_path / "extracted_domains"
    domains_path.mkdir(exist_ok=True)
    seqs_path = out_path / "sequences"
    seqs_path.mkdir(exist_ok=True)

    # HMM models
    hmm_paths = Path("../../data/raw/pfam_hmms")

    ###########
    # Action! #
    ###########

    columns = ["Uniprot ID", "PDB ID", "Chain ID", "Domain ID",
               "Domain start", "Domain end", "Coverage"]
    data = []
    prot_to_phospho_positions = {}
    parser = PDBParser(QUIET=True)

    print("Retrieving and processing domain structures...")
    for pfam_domain in tqdm(unique_pfam_domains):
        #print(f"Retrieving structures for {pfam_domain}")

        domain_pfam_mapping = pdb_pfam_mapping.loc[pdb_pfam_mapping["PFAM_ACCESSION"]==pfam_domain]
        domain_chain_pfam = pdb_chain_pfam.loc[pdb_chain_pfam["PFAM_ID"]==pfam_domain]
        
        print(f"{len(domain_pfam_mapping)} hits found for domain {pfam_domain}")
        # Create subdirectories for this particular domain
        full_structure_path = full_path / pfam_domain
        full_structure_path.mkdir(exist_ok=True)
        extracted_chain_path = chains_path / pfam_domain
        extracted_chain_path.mkdir(exist_ok=True)
        extracted_domain_path = domains_path / pfam_domain
        extracted_domain_path.mkdir(exist_ok=True)
        domain_seqs_path = seqs_path / pfam_domain
        domain_seqs_path.mkdir(exist_ok=True)

        # Some proteins have more than one copy of a domain
        # This dict keeps track of how many copies of domain
        # we have gone through, if there are multiple hits
        index_tracker = defaultdict(int)

        for idx, row in tqdm(domain_pfam_mapping.iterrows(),desc=f"{pfam_domain}"):

            pdb_code = row["PDB"]
            chain_id = row["CHAIN"]
            uniprot_id = row["UNIPROT_ACCESSION"]
            domain_id = row["PFAM_ACCESSION"]
            domain_start = row["AUTH_PDBRES_START"]
            domain_end = row["AUTH_PDBRES_END"]

            match = domain_chain_pfam.loc[(domain_chain_pfam["PDB"]==pdb_code) & (domain_chain_pfam["CHAIN"]==chain_id)]
            if len(match) > 1:
                index = index_tracker[f"{pdb_code}_{chain_id}"]
                match = match.iloc[index]
                domain_coverage = match.COVERAGE
                index_tracker[f"{pdb_code}_{chain_id}"] += 1
            else:
                domain_coverage = float(match["COVERAGE"])
                
            pdb_out_path = str(full_structure_path / f"{pdb_code}.pdb")
            chain_out_path = str(extracted_chain_path / f"{pdb_code}_{chain_id}.pdb")
            domain_out_path = str(extracted_domain_path / f"{pdb_code}_{chain_id}.pdb")

            ######################################
            # Retrieve and process the structure #
            ######################################
            if not Path(domain_out_path).is_file():
                pdb_success = query_pdb(pdb_code, full_structure_path)

                if not pdb_success:
                    warn(f"Error retrieving {pdb_code} for domain {pfam_domain}", RuntimeWarning)
                    continue

                # Extract chain and clean it up
                #print(f"Extracting chain for {pdb_code},{chain_id}...")
                call_oneinall(pdb_out_path, chain_id, chain_out_path)

                #print(f"Extracting domain {pfam_domain} from {pdb_code},{chain_id}...")
                try:
                    call_selres_range(chain_out_path, domain_out_path,
                                    domain_start, domain_end)
                except ValueError:
                    warn(f"pdb_selres failed unexpectedly on {pdb_code}_{chain_id} (range {domain_start}-{domain_end})")
                    # Delete empty file
                    domain_out_path.unlink()
                    continue
                
                # #################################
                # # Find phosphorylated positions #
                # #################################

                # print("Finding phosphosites...")
                # # TODO: ensure this is the DOMAIN only
                structure = parser.get_structure(pdb_code, domain_out_path)
                chain = structure[0][chain_id]
                # phospho_positions = find_phospho_positions(chain)
                # key = f"{pdb_code}_{chain_id}_{pfam_domain}"
                # prot_to_phospho_positions[key] = phospho_positions

                ##############################################
                # Extract sequence from the domain structure #
                ##############################################

                domain_sequence = extract_sequence_from_chain(chain)
                fasta_dict = {f"{pdb_code}_{chain_id}":domain_sequence}
                fasta_output = domain_seqs_path / f"{pdb_code}_{chain_id}.fasta"
                write_fasta(fasta_dict, fasta_output)

            if Path(domain_out_path).is_file():
                pdb_data = [uniprot_id, pdb_code, chain_id, domain_id,
                            domain_start, domain_end,
                            str(domain_coverage)]
                data.append(pdb_data)


    #Write data about the individual structures
    data = pd.DataFrame(data, columns=columns)
    data.to_csv(str(out_path / "domain_data.tsv"), sep="\t",index=None)

    # Dataframe with phosphorylation information
    # phospho_df = []
    # p_columns = ["Uniprot ID","PDB ID","Chain ID","Pfam ID","Phosphorylated","Phosphosites"]
    # for pfam_domain in tqdm(unique_pfam_domains):

    #     domain_seqs_path = seqs_path / pfam_domain
    #     hmm_model_path = hmm_paths / f"{pfam_domain}.hmm"
    #     domain_alignment = align_domain_sequences(domain_seqs_path, hmm_model_path)
    #     seq_maps = map_seqs_to_alignment(domain_alignment)
    #     # Map phosphosites to the alignment!
    #     extracted_domain_path = domains_path / pfam_domain
    #     find_phospho_idxs(extracted_domains_path, seq_maps)



    print("Et voila!")
