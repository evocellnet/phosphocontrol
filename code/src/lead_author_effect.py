"""
"""

import pandas as pd
from tqdm.notebook import tqdm

from pathlib import Path
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from tqdm import tqdm

def get_rmsd_summary(df_paths, structure_pair_to_same_author=None):

    rmsd_summary_df = []
    between_medians, within_np_medians, within_p_medians, medians_per_psite = calculate_medians(df_paths, structure_pair_to_same_author)

    columns = ["protein","residue","phosphosite","rmsd_between_groups","rmsd_within_phospho","rmsd_within_nonphospho"]
    for psite, v in medians_per_psite.items():
        protein, residue = psite.split("_")
        between_groups, within_phospho, within_nonphospho = v
        row = [protein, residue, psite, between_groups, within_phospho, within_nonphospho]
        rmsd_summary_df.append(row)

    rmsd_summary_df = pd.DataFrame(rmsd_summary_df, columns=columns)
    return rmsd_summary_df

def calculate_medians(df_paths, structure_pair_to_same_author):
    
    between_medians = []
    within_np_medians = []
    within_p_medians = []
    
    medians_per_psite = {}
    
    for rmsd_df in df_paths:
        psite = rmsd_df.parts[-2]
        df = pd.read_csv(rmsd_df,dtype={"PDB_ID_A":str,"PDB_ID_B":str,
                                        "Chain_A":str, "Chain_B":str})
        
        df["structure_a"] = df["PDB_ID_A"] + "_" + df["Chain_A"]
        df["structure_b"] = df["PDB_ID_B"] + "_" + df["Chain_B"]

        if structure_pair_to_same_author:
            drop_idxs = []
            for idx, row in df.iterrows():
                structure_a = row["structure_a"]
                structure_b = row["structure_b"]
                same_author = structure_pair_to_same_author[(structure_a, structure_b)]
                if same_author == 1:
                    drop_idxs.append(idx)
            df = df.drop(index=drop_idxs)
        
        btwn = df[df["Group"] == "between_groups"]
            
        btwn_rmsd = btwn["RMSD"].values
        between_medians.append(np.nanmedian(btwn_rmsd))

        within_np = df[df["Group"] == "within_nonphospho"]
        within_np_rmsd = within_np["RMSD"].values
        within_np_medians.append(np.nanmedian(within_np_rmsd))

        within_p = df[df["Group"] == "within_phospho"]
        within_p_rmsd = within_p["RMSD"].values
        within_p_medians.append(np.nanmedian(within_p_rmsd))
        
        medians_per_psite[psite] = (np.nanmedian(btwn_rmsd),np.nanmedian(within_p_rmsd),np.nanmedian(within_np_rmsd))
    
    return between_medians, within_np_medians, within_p_medians, medians_per_psite

def plot_rmsd_boxplot(rmsd_summary_df, out_path):
    """
    """
    sns.set_style('white')
    sns.set_context('paper',rc={"xtick.labelsize":14,"ytick.labelsize":14})
    ax = sns.boxplot(data=rmsd_summary_df, orient="h", notch=True, showfliers = False)
    ax = sns.stripplot(data=rmsd_summary_df, orient="h",palette='dark:black',alpha=0.2,jitter=True)

    ax.set_yticklabels(["Comparison between groups","Within phosphorylated proteins",
                        "Within non-phosphorylated proteins"])
    ax.set_xlabel("Median backbone RMSD (Å)",fontsize=18)
    plt.title('Excluding comparisons with the same lead author')

    plt.xlim(0,10)
    sns.despine(offset=20)
    plt.savefig(out_path / "boxplot_groups_nosameauthor.png",dpi=150,bbox_inches='tight')

def plot_rmsd_scatterplot(rmsd_summary_df, out_path):
    """
    """
    
    sns.set_context('paper',rc={"xtick.labelsize":14,"ytick.labelsize":14})
    fig = plt.figure(figsize=(8,8))
    axs = sns.scatterplot(data=rmsd_summary_df, x="rmsd_within_phospho",y="rmsd_within_nonphospho",
                    alpha=0.35,c="tab:blue",s=100)
    add_identity(axs, color='gray', ls='--', lw=1)
    plt.ylabel("Median backbone RMSD (Å) \n(non-phosphorylated structures)",fontsize=18)
    plt.xlabel("Median backbone RMSD (Å) \n(phosphorylated structures)",fontsize=18)
    plt.xlim(-.5,12.5)
    plt.ylim(-.5,12.5)

    sns.despine(offset=10)
    plt.savefig(out_path / "scatter_rmsd_phospho_vs_nonphospho_nosame.png",dpi=150, bbox_inches='tight')

def add_identity(axes, *line_args, **line_kwargs):
    identity, = axes.plot([], [], *line_args, **line_kwargs)
    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes


if __name__ == "__main__":

    #########
    # Setup #
    #########

    out_path = Path("../../results/rmsd_excluding_same_author")
    out_path.mkdir(exist_ok=True)

    glm_df = pd.read_csv("../../data/processed/glm_control_input/glm_df.csv")
    print('Building dictionary structure pair - same author..')
    structure_pair_to_same_author = {}
    for idx, row in tqdm(glm_df.iterrows()):
        structure_a, structure_b = row["structure_pair"].split("$")
        same_lead_author = row["same_lead_author"]
        structure_pair_to_same_author[(structure_a, structure_b)] = same_lead_author

    print('Excluding comparisons involving the same author...')
    rmsds_path = Path("../../results/rmsds")
    rmsd_dfs = list(rmsds_path.glob('**/rmsds_df.csv'))
    rmsd_summary_df_nosame = get_rmsd_summary(rmsd_dfs, structure_pair_to_same_author)
    rmsd_summary_df_nosame.to_csv(out_path / "rmsd_summary_nosameauthor.csv")

    # Exclude duplicates
    rmsd_summary_df_nosame_nodups = rmsd_summary_df_nosame.drop_duplicates(subset=list(rmsd_summary_df_nosame.columns)[3:])

    plot_rmsd_boxplot(rmsd_summary_df_nosame_nodups, out_path)
    plot_rmsd_scatterplot(rmsd_summary_df_nosame_nodups, out_path)
    
    print('Et voila!')