#################################################################
# Script to write eigenvector matrices for a given NMA analysis #
#################################################################
library(bio3d)

# Parse CLI arguments
args = commandArgs(trailingOnly=TRUE)
input.path <- args[1]
out.path <- args[2]
# Load modes file
load(input.path)

# Function to write each slice to a CSV file
write_slices_to_csv <- function(matrix_3d, out.path, cleaned_row_names, file_prefix = "eigenvectors") {
  dims <- dim(matrix_3d)
  O <- dims[3]
  
  for (i in 1:O) {
    slice <- matrix_3d[, , i]
    file_name <- paste0(file_prefix, "_", i, ".csv")
    
    row_name <- cleaned_row_names[i]
    file_name <- paste0(file_prefix, "_", row_name, ".csv")
    file_path <- file.path(out.path, file_name)
    write.csv(slice, file_path, row.names = FALSE)
  }
}


row_names <- rownames(modes$rmsip)
cleaned_row_names <- sub("\\.pdb$", "", row_names)

write_slices_to_csv(modes$U.subspace, out.path, cleaned_row_names)
print("Done writing matrices!")

