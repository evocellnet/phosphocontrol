# Global comparative structural analysis of responses to protein phosphorylation

Project on studying the structural effects of protein phosphorylation.
Contains code for the manuscript [Global comparative structural analysis of responses to protein phosphorylation](https://doi.org/10.1101/2024.10.18.617420).

# Software requirements

## OS requirements
The code has been run and tested on Linux (Ubuntu 20.04).

## Python dependencies
The code was run on Python v3.10.13. The code depends on the following libraries (largely the standard Python scientific stack):

- `numpy`  (used v1.23.4)
- `scipy`  (used v1.9.3)
- `numba`  (used v0.56.3)
- `pandas`  (used v1.5.2)
- `scikit-learn`  (used v1.2.0)
- `tqdm`  (used v4.64.1)
- `hdbscan`  (used v0.8.29)
- `matplotlib`  (used v3.3.4)
- `seaborn`  (used v0.12.2)
- `biopython`  (used v1.76)
- `prody`  (used v2.4.0)
- `geometricus`  (used v0.1.2)
- `PSA` (used v1.1)

## R dependencies
The code was run on R v4.4.1.

- `bio3d` (used v2.4-4)
- `FSA`  (used v0.9.5)

## External dependencies
- `Clustal Omega` (used v1.2.4)

## Installation

Clone this repository to your local phosphocontrol folder:

    cd ~ && mkdir -p phosphocontrol && cd phosphocontrol
    git clone https://github.com/evocellnet/phosphocontrol.git

# Dataset
For convenience, some data is available in the repository (paired structures dataset, functional scores, UniProt site annotation...). Due to its size, the used structure data is available instead as a [Zenodo repository](https://doi.org/10.5281/zenodo.14217157). For testing purposes, a small sample of the dataset of paired structures is provided in data/processed/pdb_pairs/demo_filtered_df.csv.

# Instructions
Individual scripts are typically set up to run on the data without the need for user input; instructions on the necessary inputs are otherwise provided as part of the documentation.

# License
The code is available under the BSD-3 license. See the LICENSE file for details.

# Citation

```bash
@article{correa2024global,
  title={Global comparative structural analysis of responses to protein phosphorylation},
  author={Correa Marrero, Miguel and Mello, Victor Hugo and Sartori, Pablo and Beltrao, Pedro},
  journal={bioRxiv},
  pages={2024--10},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

