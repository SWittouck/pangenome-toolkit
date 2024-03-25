# This script demonstrates the SCARAP sample module using four test runs on a
# dataset of L. plantarum genomes. Each run of the sample module produces a 
# genome sampling order and an ANI matrix. The four test runs use different ANI
# variants: cANI/tcANI and full core genome vs 100 core genes. cANI is the 
# normal core genome ANI, while tcANI is a trimmed version where 5% of the core
# genes are removed on both sides of the per-gene ANI distribution before cANI
# calculation. 

# dependencies: R v4.2.3, tidyverse v1.3.1, ggtree v3.8.0

library(tidyverse)
library(ggtree)
source("src/03_demonstrate_sample_functions.R")

#########################################
# Read samplings, ANI matrices and tree #
#########################################

# specify paths 
din <- "data/benchmark_scarap_v3/scarap_sample/lactiplantibacillus"
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

# explore the correlation between cANI from full core genome and subset
contrasts_canis <- 
  novcurves %>%
  filter(dataset_sampling == "cani_full", dataset_dist == "cani_full") %>%
  filter(! is.na(reference)) %>%
  mutate(
    dist_full = novelty,
    dist_sub = map2_dbl(seed, reference, ~ distmats$cani_sub[[.x, .y]])
  )
contrasts_canis %>%
  # filter(dist_full != 0, dist_sub != 0) %>%
  mutate(across(c(dist_full, dist_sub), ~ {.[. == 0] <- 1e-6 ; .})) %>%
  {cor.test(log10(.$dist_full), log10(.$dist_sub))}
(
  fig_corr_canis <-
    contrasts_canis %>%
    ggplot(aes(x = dist_full, y = dist_sub)) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    geom_point(size = 1) + 
    scale_x_log10() + scale_y_log10() + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    xlab("novelty (all SCOs)") + 
    ylab("novelty (100 SCOs)") +
    geom_text("cANI-based novelty")
)

# explore the correlation between tcANI from full core genome and subset
contrasts_tcanis <- 
  novcurves %>%
  filter(dataset_sampling == "tcani_full", dataset_dist == "tcani_full") %>%
  filter(! is.na(reference)) %>%
  mutate(
    dist_full = novelty,
    dist_sub = map2_dbl(seed, reference, ~ distmats$tcani_sub[[.x, .y]])
  )
contrasts_tcanis %>%
  # filter(dist_full != 0, dist_sub != 0) %>%
  mutate(across(c(dist_full, dist_sub), ~ {.[. == 0] <- 1e-6 ; .})) %>%
  {cor.test(log10(.$dist_full), log10(.$dist_sub))}
(
  fig_corr_tcanis <- 
    contrasts_tcanis %>%
    ggplot(aes(x = dist_full, y = dist_sub)) +
    geom_point(size = 1) +
    geom_abline(slope = 1, intercept = 0, lty = 2) + 
    scale_x_log10() + scale_y_log10() + 
    theme_bw() +
    xlab("1 - tcANI (full core genome)") + ylab("1 - tcANI (100 core genes)")
)

# visualize novelty curves (apart) 
figs_novelties <- 
  novcurves %>%
  filter(dataset_sampling == dataset_dist) %>%
  split(.$dataset_sampling) %>%
  map(
    ~ ggplot(., aes(x = step, y = novelty)) +
      geom_point(size = 0.4) + 
      scale_y_log10() +
      theme_bw()
  )

# make figure panel with novelty curves and correlation scatterplots 
fig <- function(letter, plot) {
  g <- ggplotGrob(plot + ggtitle(letter))
  g$layout$l[g$layout$name == "title"] <- 1
  g
}
ggpubr::ggarrange(
  fig("A", figs_novelties$cani_full), 
  fig("B", figs_novelties$cani_sub),
  fig("C", figs_novelties$tcani_full), 
  fig("D", figs_novelties$tcani_sub),
  fig("E", fig_corr_canis),
  fig("F", fig_corr_tcanis),
  ncol = 2, nrow = 3
)
ggsave(
  paste0(dout, "/novelties_correlations_panel.png"), units = "cm", width = 16, 
  height = 20
)

# make separate figures for cANI full and sub
figs_novelties$cani_full
ggsave(
  paste0(dout, "/novelties_cani_full.png"), units = "cm", width = 10, 
  height = 8
)
figs_novelties$cani_sub
ggsave(
  paste0(dout, "/novelties_cani_sub.png"), units = "cm", width = 10, 
  height = 8
)

# visualize novelty curves in a facet grid
corelabs <- c("full core genome", "100 core genes")
novcurves %>%
  separate(
    dataset_sampling, into = c("ani_type", "core_genes"), sep = "_", remove = F
  ) %>%
  mutate(ani_type = factor(
    ani_type, levels = c("cani", "tcani"), labels = c("cANI", "tcANI")
  )) %>%
  mutate(core_genes = factor(
    core_genes, levels = c("full", "sub"), labels = corelabs
  )) %>%
  mutate(dist = str_extract(dataset_dist, "[^_]+$")) %>%
  ggplot(aes(x = step, y = novelty, col = dist)) +
  geom_point(size = 0.4) + 
  facet_grid(rows = vars(ani_type), cols = vars(core_genes)) +
  scale_y_log10() +
  theme_bw() +
  scale_color_brewer(palette = "Paired") 
ggsave(
  paste0(dout, "/novelties_grid.png"), units = "cm", width = 16, height = 16
)

##################################
# Previous figures but as facets #
##################################

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
  paste0(dout, "/novelties_correlations_panel2.png"), units = "cm", width = 16, 
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
tree %>%
  tidygenomes::add_rootbranch() %>%
  ggtree::ggtree(layout = "rectangular") %<+%
  sampling +
  geom_tiplab2(aes(label = step), size = 1)
ggsave(
  paste0(dout, "/tree_rectangle.png"), units = "cm", width = 12, height = 50
)
