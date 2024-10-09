#!/usr/bin/env Rscript

# This script evaluates the SCARAP core module, by making a comparison with the
# full pangenome on a dataset with Lactiplantibacillus genomes.
# 
# Abbreviations: 
# - coremod = core module 
# - panmod = pan module

# dependencies: R v4.2.3, tidyverse v1.3.1, eulerr v6.1.1

library(tidyverse)
source("src/02_demonstrate_core_functions.R")

##################################
# Read core genome and pangenome #
##################################

# define relevant paths 
fin_core <- "data/benchmark_scarap_v4/scarap_core/lactiplantibacillus/corefull/genes.tsv"
fin_pan <- "data/benchmark_scarap_v4/scarap_core/lactiplantibacillus/pan/pangenome.tsv"
dout <- "results/core_demonstration"

# create output  folder 
if (! dir.exists(dout)) dir.create(dout)

# read and join datasets 
colnames <- c("gene", "genome", "orthogroup_coremod")
core <- fin_core %>% read_tsv(col_types = cols(), col_names = colnames)
colnames <- c("gene", "genome", "orthogroup_panmod")
pan <- fin_pan %>% read_tsv(col_types = cols(), col_names = colnames)

# join datasets
genes <- pan %>% left_join(core, by = c("gene", "genome"))

# remove unneeded objects
rm(colnames, core, pan, fin_core, fin_pan)

#####################################
# Compare core module to pan module #
#####################################

# calculate the number of genomes 
n_genomes <- length(unique(genes$genome))

# calculate single-copy core orthogroups (SCOs) in pan module 
scos_pan <-
  genes %>%
  group_by(orthogroup_panmod) %>%
  summarize(
    genes_panmod = n(), occurrence_panmod = sum(table(genome) == 1), 
    .groups = "drop"
  ) %>%
  filter(occurrence_panmod >= 0.95 * {{n_genomes}}) %>%
  pull(orthogroup_panmod)

# calculate SCOs in core module 
scos_core <- unique(genes$orthogroup_coremod) %>% discard(is.na)

# determine mapping between pan and core SCOs 
scos_both <- 
  genes %>%
  filter(orthogroup_panmod %in% scos_pan, orthogroup_coremod %in% scos_core) %>%
  group_by(orthogroup_panmod, orthogroup_coremod) %>%
  summarize(n_genomes = length(unique(genome)), .groups = "drop") %>%
  filter(n_genomes >= 0.95 * {{n_genomes}})

# calculate some statistics 
numbers <-
  list(
    scos_core = length(scos_core),
    scos_pan = length(scos_pan), 
    scos_both = nrow(scos_both), 
    precision = nrow(scos_both) / length(scos_core),
    recall = nrow(scos_both) / length(scos_pan)
  )
numbers$f_measure = 2 / (1 / numbers$precision + 1 / numbers$recall)

# write statistics
capture.output(numbers, file = paste0(dout, "/statistics.txt"))

# all pan and core SCOs should be in the shared SCOs
all(scos_both$orthogroup_panmod %in% scos_pan)
all(scos_both$orthogroup_coremod %in% scos_core)
