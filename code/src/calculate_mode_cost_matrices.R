library(bio3d)
library(lsa)

#######################
# Parse CLI arguments #
#######################
args = commandArgs(trailingOnly=TRUE)

#########
# Setup #
#########

input.path <- args[1]
out.path <- args[2]
# Load modes file
load(input.path)

m <- dim(modes$U.subspace)[1]
n <- dim(modes$U.subspace)[2]
o <- dim(modes$U.subspace)[3]

# Compute all unique pairs of elements in the third dimension
combinations <- combn(o, 2)

# Initialize a list to store the cost matrices
#cost_matrices <- array(0, dim = c(n, n, ncol(combinations)))

# Compute pairwise similarities for each combination of the third dimension
for (k in 1:ncol(combinations)) {
  a <- combinations[1, k]
  b <- combinations[2, k]
  #cat("Computing similarity matrix for combination (", a, ", ", b, ")...\n", sep = "")
  # Compute the similarity matrix for the combination (a, b)
  similarity_matrix <- outer(
    1:n, 1:n,
    Vectorize(function(i, j) cosine(modes$U.subspace[, i, a], modes$U.subspace[, j, b]))
  )
  cost_matrix <- 1 - similarity_matrix
  #cost_matrices[, , k] <- cost_matrix
  
  # Write the similarity matrix to a CSV file in the output directory
  # The indexes in the filenames are 0-indexed
  file_name <- file.path(out.path, paste("cost_matrix_", a-1, "_", b-1, ".csv", sep = ""))
  # cat("Writing similarity matrix to file: ", file_name, "\n", sep = "")
  write.csv(cost_matrix, file = file_name, row.names = FALSE)
}


# Initialize the cost matrix
#cost_matrix <- matrix(0, n, n)

# Compute cost matrices and solve the assignment problem for each pair of samples
#for (a in 1:(p-1)) {
#  for (b in (a+1):p) {
#    # Initialize the cost matrix for the pair (a, b)
#    cost_matrix <- matrix(0, n, n)
    
#    # Compute the cost matrix
#    for (i in 1:n) {
#      for (j in 1:n) {
#        cost_matrix[i, j] <- 1 - cosine(modes$U.subspace[, i, a], modes$U.subspace[, j, b])
#      }
#    }
#    filename <- paste0("cost_matrix_", a-1, "_", b-1, ".csv")
#    file_path <- file.path(out.path, filename)
#    write.csv(cost_matrix, file = file_path, row.names = FALSE)
#  }
#}    
