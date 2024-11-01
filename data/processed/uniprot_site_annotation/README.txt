uniprot_site_annotation.csv

Each row corresponds to an annotated residue. 

uniprot_id:   UniProt ID
site_type:    either Active site, DNA binding site, Binding site or Site.
              Active site:  residues directly involved in catalysis.
                            If the role of the residue is known, this is indicated
                            using a controlled vocabullary (see 'description')
              Binding site: Binding site for any chemical group (metals, cofactors,
                            substrates and products of enzymes, ligands...)
                            The chemical group is described in 'ligand'
              DNA binding:  Residues annotated to bind DNA.
              Site:         Interesting residues not defined in other sections.
                            Role described in 'description'.
residue:      residue index (according to the UniProt sequence)
description:  role performed by Active site or Site residues
ligand:       chemical group binding to Binding site residues
