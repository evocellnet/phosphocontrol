# Script to test for differences in changes in RMSD upon phosphorylation by residue 

library(FSA)
library(here)

data_path <- here("notebooks","backbone_comparison","rmsd_merged_df.csv")
df <- read.csv(data_path)

kruskal.test(rmsd_between_groups~CHEM_COMP_ONE, data= df)

dt = dunnTest(rmsd_between_groups ~ CHEM_COMP_ONE,
              data=df,
              method="bh")
