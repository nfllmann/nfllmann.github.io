## Functions for set operations ##

# Calculate the symmetric difference
vol_diff_sym_x_fixed <- function(ij, is_in_gamma_u) {
  # Computes the symmetric difference of Gamma_i and Gamma_j 
  # Parameters:
  #   ij: Indices of the elements to compute the symmetric difference
  #   is_in_gamma_u: Matrix representing random sets
  # Returns:
  #   Symmetric difference
  return(mean(is_in_gamma_u[,ij[1]] | is_in_gamma_u[,ij[2]]) - mean(is_in_gamma_u[,ij[1]] & is_in_gamma_u[,ij[2]]))
}

# Parallel computation of k_{set}
rbf_set_hsic_parall_x_fixed <- function(is_in_gamma_u, n_u, param = NULL) {
  # Computes the Gram matrix of a random set, i.e., the matrix of coefficients K(Gamma_i, Gamma_j)
  # Parameters:
  #   is_in_gamma_u: Matrix representing a random set
  #   n_u: Number of set realizations
  #   param: Optional parameter for kernel calculation
  # Returns:
  #   Gram matrix 
  
  # Initialize the Gram matrix
  gram_matrix <- matrix(0, n_u, n_u)
  
  # Generate pairs of indices 
  indice <- matrix(0, ncol = 2, nrow = (n_u - 1) * n_u / 2)
  k <- 0
  for (i in 1:(n_u - 1)) {
    for (j in (i + 1):n_u) {
      k <- k + 1
      indice[k, ] <- c(i, j)
    }
  }
  
  # Parallel computation of the symmetric volume difference
  n.cores <- detectCores()
  clust <- makeCluster(n.cores)
  clusterExport(clust, c("vol_diff_sym_x_fixed"))
  tempo <- parApply(clust, indice, 1, vol_diff_sym_x_fixed, is_in_gamma_u = is_in_gamma_u)
  stopCluster(clust)  # Stop cluster after computation
  
  # Fill the Gram matrix with computed values
  k <- 0
  for (i in 1:(n_u - 1)) {
    for (j in (i + 1):n_u) {
      k <- k + 1
      gram_matrix[i, j] <- tempo[k]
      gram_matrix[j, i] <- gram_matrix[i, j]
    }
  }
  
  # Compute parameter if not provided
  if (is.null(param)) {
    param <- sqrt(mean(gram_matrix[upper.tri(gram_matrix)]))
  }
  
  # Apply the RBF kernel function to the Gram matrix
  gram_matrix <- exp(-gram_matrix / (2 * param^2))
  
  return(gram_matrix)
}
