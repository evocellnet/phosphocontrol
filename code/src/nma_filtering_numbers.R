#!/usr/bin/env Rscript
# Set number of cores to use
library(bio3d)
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

phospho_nr_before_filtering <- sum(annot$Status == "Phosphorylated")
nonphospho_nr_before_filtering <- sum(annot$Status == "Non-phosphorylated")

# Remove structures with missing residues
omitted_nr <- length(which(!conn))
total_nr <- length(which(conn)) + length(which(!conn))
print(paste0("Total number of structures: ",total_nr))
print(paste0("Structures with missing residues: ",omitted_nr))

# Write to file number of starting structures
# Number of phosphorylated structures left after filtering
# Number of non-phosphorylated structures left after filtering

phospho_nr_after_filtering <- 0
nonphospho_nr_after_filtering <- 0

try({
  pdbs <- trim(pdbs, row.inds=which(conn))
  
  excluded_names <- names(which(!conn))
  excluded_names <- gsub(".*/","",excluded_names)
  excluded_names <- gsub("\\..*","", excluded_names)
  annot <- annot[!(annot$Full.ID %in% excluded_names), , drop = FALSE]
  
  phospho_nr_after_filtering <- length(which(annot$Status == "Phosphorylated"))
  nonphospho_nr_after_filtering <- length(which(annot$Status == "Non-phosphorylated"))
}, silent = TRUE)

# Check if an error was raised and values are still 0
if (phospho_nr_after_filtering == 0 & nonphospho_nr_after_filtering == 0) {
  message("There are no phosphorylated and non-phosphorylated structures left after trimming")
}

# pdbs <- trim(pdbs, row.inds=which(conn))

# excluded_names <- names(which(!conn))
# excluded_names <- gsub(".*/","",excluded_names)
# excluded_names <- gsub("\\..*","", excluded_names)
# annot <- annot[!(annot$Full.ID %in% excluded_names), , drop = FALSE]

# # Exit if we now have no phosphorylated or non-phosphorylated structures
# phospho_nr_after_filtering <- length(which(annot$Status == "Phosphorylated"))
# nonphospho_nr_after_filtering <- length(which(annot$Status == "Non-phosphorylated"))

df <- data.frame(Number = c(total_nr, phospho_nr_before_filtering, nonphospho_nr_before_filtering, phospho_nr_after_filtering, nonphospho_nr_after_filtering))
rownames(df) <- c("total_number", "phospho_before", "nonphospho_before", "phospho_after", "nonphospho_after")
out_df_path <- file.path(out.path, "structure_numbers.csv")
write.csv(df, file = out_df_path, row.names = TRUE)

