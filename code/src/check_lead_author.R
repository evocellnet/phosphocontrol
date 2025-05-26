
# This script fits generalized linear models (GLMs) to investigate how phosphorylation status
# and lead authorship influence structural similarity between protein pairs, as measured by RMSD.
# Three different GLMs are fit:
#   1. Across all protein pairs (including both phosphorylated and non-phosphorylated)
#   2. Subset to only phosphorylated protein pairs
#   3. Subset to only non-phosphorylated protein pairs
#
# Usage: this script is meant to be used interactively.
#
# Inputs:
#   - "../../data/processed/glm_control_input/glm_df.csv" 
#       Contains all protein pairs with columns:
#         - rmsd: numeric, root mean square deviation between structures
#         - is_phospho_a: factor, whether protein A is phosphorylated
#         - is_phospho_b: factor, whether protein B is phosphorylated
#         - same_lead_author: factor, whether both structures are from the same lead author
#
#   - "../../data/processed/glm_control_input/phospho_df.csv" 
#       Contains only phosphorylated protein pairs (subset of the above)
#
#   - "../../data/processed/glm_control_input/nonphospho_df.csv" 
#       Contains only non-phosphorylated protein pairs (subset of the above)
#
# Output: GLM model objects (model, p_model, np_model) held in memory for further inspection.
#         Models may be inspected using `summary(model)`.
#
# Note: the working directory must be set correctly for relative file paths to resolve.

# Set working directory to script location
setwd("~/Documents/projects/phosphocontrol/code/src")

# All protein pairs
df <- read.table("../../data/processed/glm_control_input/glm_df.csv",sep=",",header=TRUE)

df$is_phospho_a <- as.factor(df$is_phospho_a)
df$is_phospho_b <- as.factor(df$is_phospho_b)
df$same_lead_author <- as.factor(df$same_lead_author)

model <- glm(rmsd ~ is_phospho_a + is_phospho_b + same_lead_author, family="gaussian", data=df)

# Phosphorylated protein pairs
p_df <- read.table("../../data/processed/glm_control_input/phospho_df.csv",sep=",",header=TRUE)
p_df$is_phospho_a <- as.factor(p_df$is_phospho_a)
p_df$is_phospho_b <- as.factor(p_df$is_phospho_b)
p_df$same_lead_author <- as.factor(p_df$same_lead_author)

p_model <- glm(rmsd ~ same_lead_author, family="gaussian", data=p_df)

# Non-phosphorylated protein pairs
np_df <- read.table("../../data/processed/glm_control_input/nonphospho_df.csv",sep=",",header=TRUE)
np_df$is_phospho_a <- as.factor(np_df$is_phospho_a)
np_df$is_phospho_b <- as.factor(np_df$is_phospho_b)
np_df$same_lead_author <- as.factor(np_df$same_lead_author)

np_model <- glm(rmsd ~ same_lead_author, family="gaussian",data=np_df)

