#!/usr/bin/env Rscript

# This script demonstrates the SCARAP sample module using four test runs on a
# dataset of Lactiplantibacillus genomes. Each run of the sample module produces a 
# genome sampling order and an ANI matrix. The four test runs use different ANI
# variants: cANI/tcANI and full core genome vs 100 core genes. cANI is the 
# normal core genome ANI, while tcANI is a trimmed version where 5% of the core
# genes are removed on both sides of the per-gene ANI distribution before cANI
# calculation. 

# dependencies: R v4.2.3, tidyverse v1.3.1, ggtree v3.8.0, ape v5.8

library(tidyverse)
library(ggtree)
source("src/03_demonstrate_sample_functions.R")

#########################################
# Read samplings, ANI matrices and tree #
#########################################

# specify paths 
din <- "data/benchmark_scarap_v4/scarap_sample/lactiplantibacillus"
dout <- "results/sample_demonstration"

# make output folder 
if (! dir.exists(dout)) dir.create(dout)

# define the names and folder names of the datasets
datasets <- 
  c("corefull_mean", "core100_mean", "corefull_mean90", "core100_mean90") %>%
  set_names(c("cani_full", "cani_sub", "tcani_full", "tcani_sub"))

# read samplings for the four datasets
samplings <- 
  datasets %>%
  map_chr(~ paste0(din, "/", ., "/seeds.txt")) %>%
  map(read_lines)

# read ANI tables for the four datasets and convert to distance matrices
distmats <- 
  datasets %>%
  map_chr(~ paste0(din, "/", ., "/identities.tsv")) %>%
  map(read_tsv, col_types = cols()) %>%
  map(tibble2matrix) %>%
  map(~ 1 - .)

# read tree
tree <- paste0(din, "/tree/lactiplantibacillus.treefile") %>% ape::read.tree()

############################
# Calculate novelty curves #
############################

# calculate novelty curves
# dataset_sampling = which genome sampling order do we use? 
# dataset_dist = which dataset do we use to calculate the novelties? 
novcurves <- 
  tribble(
    ~ dataset_sampling, ~ dataset_dist, 
    "cani_full", "cani_full",
    "cani_sub", "cani_sub",
    "cani_sub", "cani_full",
    "tcani_full", "tcani_full",
    "tcani_sub", "tcani_sub",
    "tcani_sub", "tcani_full"
  ) %>%
  mutate(novelties = map2(
    dataset_sampling, dataset_dist, 
    ~ determine_novelties(samplings[[.x]], distmats[[.y]])
  )) %>%
  unnest(novelties)

#################################################
# Visualize novelty curves and ANI correlations #
#################################################

# extract canis for independent contrasts
contrasts_canis <- 
  novcurves %>%
  filter(dataset_sampling == "cani_full", dataset_dist == "cani_full") %>%
  filter(! is.na(reference)) %>%
  mutate(
    dist_full = novelty,
    dist_sub = map2_dbl(seed, reference, ~ distmats$cani_sub[[.x, .y]])
  )

# extract tcanis for independent contrasts
contrasts_tcanis <- 
  novcurves %>%
  filter(dataset_sampling == "tcani_full", dataset_dist == "tcani_full") %>%
  filter(! is.na(reference)) %>%
  mutate(
    dist_full = novelty,
    dist_sub = map2_dbl(seed, reference, ~ distmats$tcani_sub[[.x, .y]])
  )

# visualize novelty curves (facet grid)
(
  fig_novelties <- 
    novcurves %>% 
    filter(dataset_sampling == dataset_dist) %>%
    separate(dataset_dist, into = c("ani", "scos"), sep = "_") %>%
    mutate(scos = recode(scos, "full" = "all SCOs", sub = "100 SCOs")) %>% 
    mutate(scos = factor(scos, levels = c("all SCOs", "100 SCOs"))) %>%
    mutate(ani = recode(ani, "cani" = "cANI", "tcani" = "tcANI")) %>%
    ggplot(., aes(x = step, y = novelty)) +
    geom_point(size = 0.4) +
    facet_grid(rows = vars(ani), cols = vars(scos)) + 
    scale_y_log10() +
    theme_bw() +
    xlab("sampling step") + ylab("novelty (1 - ANI)")
)

# explore correlations (facet wrap)
(
  fig_corr <-
    bind_rows(contrasts_canis, contrasts_tcanis) %>%
    mutate(trim = recode(
      dataset_dist, "cani_full" = "no gene trimming (cANI)",
      "tcani_full" = "gene trimming (tcANI)"
    )) %>%
    mutate(trim = fct_rev(trim)) %>%
    ggplot(aes(x = dist_full, y = dist_sub)) +
    geom_point(size = 1) +
    geom_abline(slope = 1, intercept = 0, lty = 2) + 
    facet_wrap(facets = vars(trim)) + 
    scale_x_log10() + scale_y_log10() + 
    theme_bw() +
    xlab("novelty (all SCOs)") + ylab("novelty (100 SCOs)")
)

# make figure panel with faceted figures
ggpubr::ggarrange(
  fig("A", fig_novelties), 
  fig("B", fig_corr),
  ncol = 1, nrow = 2, heights = c(1, 0.7)
)
ggsave(
  paste0(dout, "/novelties_correlations.png"), units = "cm", width = 16, 
  height = 20
)

###################################
# Visualize sampling next to tree #
###################################

# midpoint root the tree 
tree <- phytools::midpoint.root(tree)

# create table with sampling order 
sampling <- 
  tibble(genome = samplings$tcani_sub, step = 1:length(genome)) %>%
  slice(1:10)

# visualize the tree (circular)
tree %>%
  tidygenomes::add_rootbranch() %>%
  ggtree::ggtree(layout = "circular") %<+%
  sampling +
  geom_tiplab2(aes(label = step))
ggsave(paste0(dout, "/tree.png"), units = "cm", width = 12, height = 12)

# visualize the tree (rectangular)
tree2 <- tree %>% tidygenomes::add_rootbranch() 
tree2 <- 
  tree2 %>%
  ggtree::ggtree(layout = "rectangular") %<+%
  sampling +
  geom_tiplab2(aes(label = step), size = 1)
for (clade in to_collapse) {
  tree2 <- collapse(tree2, clade, "mixed")
}
tree2
ggsave(
  paste0(dout, "/tree_rectangle.png"), units = "cm", width = 12, height = 50
)

# collapse clades within the species level
to_collapse <- young_clades(tree, 0.05)
# tree2 <- tree %>% tidygenomes::add_rootbranch() 
treefig <- ggtree::ggtree(tree, layout = "rectangular")
for (clade in to_collapse) {
  treefig <- collapse(treefig, clade, "mixed")
}
treefig
ggsave(
  paste0(dout, "/tree_collapsed.png"), units = "cm", width = 12, height = 50
)
