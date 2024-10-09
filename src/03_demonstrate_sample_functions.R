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

# helper function for young_clades
check_nodes <- function(
  tree, distmat, threshold, node, to_collapse = integer()
) {
  children <- tree$edge[tree$edge[, 1] == node, 2]
  if (length(children) == 0) return(to_collapse)
  desc_left <- 
    phytools::getDescendants(tree, children[1]) %>%
    {tree$tip.label[.]} %>% 
    discard(is.na)
  desc_right <- 
    phytools::getDescendants(tree, children[2]) %>%
    {tree$tip.label[.]} %>% 
    discard(is.na)
  if (any(distmat[desc_left, desc_right] <= threshold)) {
    to_collapse <- c(to_collapse, node)
    return(to_collapse)
  }
  for (child in children) {
    to_collapse <- check_nodes(tree, distmat, threshold, child, to_collapse)
  }
  to_collapse
}

#' Determine young clades
#' 
#' This function will determine all clades in a tree that are younger than a
#' threshold value. A clade is younger than the threshold if any of the "left
#' descendants" of its LCA are closer than the threshold to any of the "right
#' descendants" of its LCA. Clades nested within young clades are not returned
#' separately. 
#' 
#' @param tree A rooted tree of class phylo. 
#' @param dists A threshold value (in terms of branch lengths). 
#'   
#' @return A vector with node names of LCAs of young clades
young_clades <- function(tree, threshold) {

  # calculate distance matrix between all tips
  distmat <- ape::cophenetic.phylo(tree)
  # assign rownames to distmat
  rownames(distmat) <- colnames(distmat)
  # determine root node 
  root <- tree$edge[, 1][! tree$edge[, 1] %in% tree$edge[, 2]][1]
  # check all nodes from root to tip
  check_nodes(tree, distmat, 0.05, root) 

}

tibble2matrix <- function(t) {
  t %>%
    as.data.frame() %>%
    `rownames<-`(.$genome) %>%
    {.[, -1]} %>%
    as.matrix()
}

fig <- function(letter, plot) {
  g <- ggplotGrob(plot + ggtitle(letter))
  g$layout$l[g$layout$name == "title"] <- 1
  g
}