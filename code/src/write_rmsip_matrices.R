########################################################################
# Script to write all pairwise RMSIP matrices for a given NMA analysis #
########################################################################
library(bio3d)

# Parse CLI arguments
args = commandArgs(trailingOnly=TRUE)
input.path <- args[1]
out.path <- args[2]
# Load modes file
load(input.path)

perform_pairwise_comparisons <- function(U_subspace, cleaned_row_names, out_path) {
  # Check dimensions
  dims <- dim(U_subspace)
  m <- dims[2] 
  n <- dims[3] # Number of elements in the third dimension
  
  # Generate all unique pairs of indices
  pairs <- combn(n, 2)
  
  # Perform pairwise comparisons
  apply(pairs, 2, function(pair) {
    i <- pair[1]
    j <- pair[2]
    
    result <- rmsip(U_subspace[,,i], U_subspace[,,j],subset=m)
    overlap <- result$overlap
        file_name <- paste0("overlap_", cleaned_row_names[i], "_vs_", cleaned_row_names[j], ".csv")
    file_path <- file.path(out_path, file_name)
    write.csv(overlap, file_path, row.names = FALSE)
  })
}

row_names <- rownames(modes$rmsip)
cleaned_row_names <- sub("\\.pdb$", "", row_names)

perform_pairwise_comparisons(modes$U.subspace, cleaned_row_names, out.path)
print("Done writing RMSIP matrices!")

