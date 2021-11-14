#' Do genotypes cluster by microhabitat?
#'
#' Calculates how many pairs of plants of the same genotype are also in the
#' same microhabitat.
#' 
#' @param genotype Vector of genotype labels.
#' @param habitat Vector of habitat labels of the same length as `genotype`.
#' 
#' @return Integer number of pairs of individuals sharing both genotype and 
#' habitat
#' @author Tom Ellis
habitat_clustering <- function(genotype, habitat){
  # if(length(genotype) != length(habitat)){
  #   stop("`genotype` and `habitat` are different lengths.")
  # }
  # For all pairs of plants, check whether genotypes are the same
  genotype_match <- combn(genotype, 2)
  genotype_match <- genotype_match[1,] == genotype_match[2,]
  # For all pairs of plants, check whether habitats are the same
  habitat_match <- combn(habitat, 2)
  habitat_match <- habitat_match[1,] == habitat_match[2,]
  # Check whether genotype and habitat match.
  both_match <- genotype_match * habitat_match
  
  return(sum(both_match, na.rm = T))
}

#' Rotate elements in a vector
#' 
#' Shift the positions of elements in vector x by n positions. Move the n
#' elements at the end of the vector to the beginning.
#' 
#' @param x Vector.
#' @param n Integer number of positions to shift the vector.
#' 
#' @return Permutation of x, shifted by n positions.
#' @author Tom Ellis
#' 
#' @example
#' rotate_vector(1:5, 3)
#' [1] 3 4 5 1 2
rotate_vector <- function(x, n){
  c(tail(x, n), head(x, -n))
}

#' Permute habitat labels
#' 
#' Reorders microhabitats, and moves them all a random amount left or right.
#' It then calculates how often plants of the same DGG are found in the same
#' microhabitat in the observed and permuted datasets.
#' 
#' @param genotype Vector of genotype labels.
#' @param habitat Vector of habitat labels of the same length as `genotype`.
#' @param nreps Number of permutations to perform
#' 
#' @return Data frame showing how often plants of the same DGG are found in the
#' same microhabitat in the observed and permuted datasets.
habitat_permutations <- function(genotype, habitat, nreps){
  obs <- habitat_clustering(genotype, habitat)
  
  permuted <- numeric(length = nreps)
  for(i in 1:nreps){
    # Split habitat labels into a list
    permute_habitat <- split(habitat, habitat)
    # Reorder habitat labels
    ix <- sample(1:length(permute_habitat), size = length(permute_habitat), replace = FALSE)
    permute_habitat <- do.call('c', permute_habitat[ix])
    # Rotate the vector by some random amount
    permute_habitat <- rotate_vector(
      x = permute_habitat,
      n = sample(1:length(permute_habitat), 1)
    )
    
    permuted[i] <- habitat_clustering(genotype, permute_habitat)
  }
  
  data.frame(
    type = c(rep("permuted", length(permuted)), "observed"),
    matches = c(permuted, obs)
  )
}