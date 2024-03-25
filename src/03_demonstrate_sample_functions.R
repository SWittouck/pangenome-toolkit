#' Determine novelty curve 
#' 
#' This function will determine the novelties of a set of sampled objects given
#' a matrix of all pairwise distances between the objects. 
#' 
#' @param seeds A list of the seeds that were sampled. 
#' @param dists A matrix with pairwise distances between the objects, with row
#'   and column names. 
#'   
#' @return A tibble with the columns step, seed and novelty
determine_novelties <- function(seeds, dists) {
  
  n <- length(seeds)
  novelties <- double(n)
  novelties[1] <- Inf
  references <- character(n)
  references[1] <- NA
  dists <- dists[seeds, seeds]
  
  for (i in 2:n) {
    dists_to_seeds <- dists[i, 1:(i - 1), drop = F] # --> one-row matrix
    min_index <- which.min(dists_to_seeds) # col number of minimum
    novelties[i] <- dists_to_seeds[1, min_index]
    references[i] <- colnames(dists_to_seeds)[min_index]
  }
  
  tibble(step = 1:n, seed = seeds, reference = references, novelty = novelties)
  
}

tibble2matrix <- function(t) {
  t %>%
    as.data.frame() %>%
    `rownames<-`(.$genome) %>%
    {.[, -1]} %>%
    as.matrix()
}
