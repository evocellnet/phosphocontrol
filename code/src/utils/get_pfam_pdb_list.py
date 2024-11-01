"""
"""

from pathlib import Path

if __name__ == "__main__":

    path = Path("../../../data/raw/pfam_structures/full_structures")
    pfam_paths = list(path.glob("**"))
    pfam_paths = [x for x in pfam_paths if x.is_dir()
                  and x.stem.startswith("PF")]

    pdb_codes = []
    for pfam_path in pfam_paths:
        pdbs = list(pfam_path.glob("*.pdb"))
        for pdb in pdbs:
            pdb_code = pdb.stem
            pdb_codes.append(pdb_code)

    print(f"Got {len(pdb_codes)} PDB codes from {len(pfam_paths)} Pfam domains")
    out_file = str(path / "all_pdb_codes.txt")
    with open(out_file, "w") as f:
        for pdb_code in pdb_codes:
            row = "".join([pdb_code, "\n"])
            f.write(row)

    print("All done here!")
