#!/usr/bin/env Rscript

# This script performs ensemble normal mode analysis (eNMA) to investigate the
# effects of phosphorylation on protein dynamics.
# It first processes a collection of PDB files, filters out structures that have
# e.g. breaks in connectivity and conducts principal component analysis of the
# structures. It then computes ensemble anisotropic network models and writes to
# disk several properties derived from this (e.g. residue fluctuations and 
# linear mutual information between residues).

# Usage:
# This script is meant to be called by dynamics_analysis_per_psite.py.
# The analysis is performed on a given folder containing structures pertaining a
# specific phosphosite.

#   Rscript nma_analysis.R $INPUT_PATH $ANNOTATION_FILE $OUTPUT_PATH
#
# Arguments
# ---------
#   $INPUT_PATH      : Directory containing PDB files to analyze for a given phosphosite.
#   $ANNOTATION_FILE : CSV file with columns 'Full.ID' (matching PDB filenames
#                      without extension)
#                      and 'Status', indicating either 'Phosphorylated' or 
#                      'Non-phosphorylated'.
#   $OUTPUT_PATH     : Directory where output files (plots, data) will be saved.
#
# Output
# ------
#   - PCA plot of atomic coordinates saved as 'pca.png'.
#   - ANM mode fluctuations plot saved as 'modes_diffs.png'.
#   - Fluctuation data saved as 'fluctuation_data.csv'.
#   - ANM modes object saved as 'modes.RData.gz'.
#   - Linear mutual information matrices saved as CSV files.
#
# Author: Miguel Correa Marrero


library(bio3d)
# Set number of cores to use
# Using more than one core leads to weird bugs sometimes
Sys.setenv(MC_CORES = 1)
print('Starting!')

#######################
# Parse CLI arguments #
#######################
args = commandArgs(trailingOnly=TRUE)

#########
# Setup #
#########

input.path <- args[1]
annot.file <- args[2]
out.path <- args[3]

dir.create(out.path)

annot <- read.table(annot.file, sep=",", header=TRUE)
# Code to remove residues with missing alpha carbon coordinates
pdb_files <- list.files(input.path, pattern="*.pdb", full.names=TRUE)

# Remove files that are not used in this comparison
annot.ids <- annot[,"Full.ID"]
# Extract filenames without path and extension from the list of all files
file_names <- gsub(".*/", "", pdb_files)
file_names <- gsub("\\..*", "", file_names)
pdb_files <- pdb_files[file_names %in% annot.ids]

for(i in pdb_files) {
  dirn <- dirname(i)
  backup <- file.path(dirn, "backup_pdb")
  pdb <- read.pdb(i)
  inds <- atom.select(pdb, "calpha")
  pat.all <- paste(pdb$atom$resid, pdb$atom$resno, pdb$atom$chain, pdb$atom$insert, sep="_")
  pat <- pat.all[inds$atom]
  inds.all <- inds
  inds.all$atom <- which(pat.all %in% pat)
  inds.all$xyz <- atom2xyz(inds.all$atom)
  pdb <- trim(pdb, inds=inds.all)
  if(!identical(pat, unique(pat.all))) {
    warning(paste("Removing residues without CA coordinates...", i))
    if(!dir.exists(backup)) {
      dir.create(backup)
    }
    file.copy(i, backup)
    write.pdb(pdb, file=i, chainter=TRUE)
  }
}

pdbs <- pdbaln(pdb_files, fit=TRUE)
pdbs <- read.all(pdbs)

conn <- inspect.connectivity(pdbs)
# Remove structures with missing residues
omitted_nr <- length(which(!conn))
total_nr <- length(which(conn)) + length(which(!conn))
print(paste0("Total number of structures: ",total_nr))
print(paste0("Structures with missing residues: ",omitted_nr))
if (omitted_nr == total_nr) {
  stop("No complete structures; exiting")
}
pdbs <- trim(pdbs, row.inds=which(conn))

excluded_names <- names(which(!conn))
excluded_names <- gsub(".*/","",excluded_names)
excluded_names <- gsub("\\..*","", excluded_names)
annot <- annot[!(annot$Full.ID %in% excluded_names), , drop = FALSE]

# Exit if we now have no phosphorylated or non-phosphorylated structures
phospho_nr <- which(annot$Status == "Phosphorylated")
nonphospho_nr <- which(annot$Status == "Non-phosphorylated")

if (length(phospho_nr)==0) {
  stop("No phosphorylated structures after filtering; exiting")
}
if (length(nonphospho_nr)==0) {
  stop("No non-phosphorylated structures after filtering; exiting")
}

# Remove conformationally redundant structures to save memory
names_phospho <- annot$Full.ID[phospho_nr]
names_nonphospho <- annot$Full.ID[nonphospho_nr]
# Separate into phosphorylated and non-phosphorylated structures
excl_np <- unlist(lapply(c(names_nonphospho), grep, pdbs$id))
excl_p <- unlist(lapply(c(names_phospho), grep, pdbs$id))
pdbs_p <- trim(pdbs, row.inds=which(!(1:length(pdbs$id) %in% excl_np)))
pdbs_np <- trim(pdbs, row.inds=which(!(1:length(pdbs$id) %in% excl_p)))

if (length(pdbs_np$id)>1) {
  redundant_np <- filter.rmsd(pdbs_np$xyz, cutoff=0.2, fit=TRUE)
  pdbs_np <- trim(pdbs_np, row.inds=redundant_np$ind)
}
if (length(pdbs_p$id)>1) {
  redundant_p <- filter.rmsd(pdbs_p$xyz, cutoff=0.2, fit=TRUE)
  pdbs_p <- trim(pdbs_p, row.inds=redundant_p$ind)
}

ids_np <- unlist(strsplit(basename(pdbs_np$id), split=".pdb"))
ids_p <- unlist(strsplit(basename(pdbs_p$id), split=".pdb"))
all_filt_ids <- c(ids_np, ids_p)
all_ids <- unlist(strsplit(basename(pdbs$id), split=".pdb"))
set_diff <- setdiff(all_ids,all_filt_ids)

pdbs <- trim(pdbs, row.inds=which(!(1:length(pdbs$id) %in% set_diff)))

#######
# PCA #
#######

# Fit PCA of atomic coordinates
print("Creating PCA of atomic coordinates...")
pc.xray <- pca(pdbs, core.find=TRUE)

pca_plot_path <- file.path(out.path, "pca.png")
png(file=pca_plot_path)
plot(pc.xray,  col=annot[,"color"])
dev.off()

###########
# Fit ANM #
###########

print("Performing elastic network calculations...")
modes <- aanma(pdbs, rtb=TRUE, reduced=TRUE,ncore=1)
modes.path <- file.path(out.path, "modes.RData.gz")
save(modes, file = modes.path, compress = "gzip")

# Save residue-wise fluctuations to a csv file
print("Saving fluctuations to disk...")
fluctuations.path <- file.path(out.path, "fluctuation_data.csv")
write.csv(modes$fluctuations, file = fluctuations.path, row.names = TRUE)

print("Plotting fluctuations...")
modes_plot_path <- file.path(out.path, "modes_diffs.png")
png(file=modes_plot_path)
plot(modes, col=annot[,"color"], signif=TRUE)
dev.off()

##########################################
# Linear mutual information calculations #
##########################################

# Calculate residue-residue linear mutual information for each structure 
# and write it to a csv
cov_matrix <- cov.enma(modes)

depth <- dim(cov_matrix)[3]
for (i in 1:depth) {
  # Access the i-th layer of the array
  cov_layer <- cov_matrix[,,i]
  # Calculate LMI for the structure
  lmi_matrix <- bio3d:::.cov2dccm(cov_layer, method="lmi")
  
  # Get the name of the structure
  split_string <- strsplit(pdb_files[i], "/")
  length_string <- length(split_string[[1]])
  pdb_name <- split_string[[1]][length_string]
  # Remove the file extension
  pdb_name <- strsplit(pdb_name, "\\.")[[1]][1]
  
  # Check if the structure is phosphorylated
  if (annot[,"color"][i]=="orange") {
    is_p <- "phospho"
  } else {
    is_p <- "nonphospho"
  }
  
  # Create filename
  file_name <- paste0("lmi_", pdb_name, "_", is_p, ".csv")
  file_path <- file.path(out.path, file_name)
  # Write csv
  write.csv(lmi_matrix, file = file_path, row.names = TRUE)
}

############################
# Similarities in dynamics #
############################

# Calculate similarities in normal modes between structures using
# the Bhattacharyya coefficient

print("Calculating Bhattacharyya coefficients...")
covs <- cov.enma(modes)
bc <- bhattacharyya(modes, covs=covs)

bcs.path <- file.path(out.path, "bcs.csv")
write.csv(bc, file = bcs.path, row.names = TRUE)

print("Et voila!")
