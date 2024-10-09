#!/usr/bin/env Rscript

# This script explores the full pangenome of a dereplicated set of > 6,000
# genomes of Lactobacillales, computed by the LEGEN pipeline.

# dependencies: R v4.2.3, tidyverse v1.3.1

library(tidyverse)

#############
# READ DATA # 
#############

# define input/output folders
din_legen <- "data/legen_v4/"
dout <- "results/legen_exploration"

# define input paths of files
fin_genomes <-
  paste0(din_legen, "/all/genomes_metadata.csv")
fin_orthogroup_occurrences <- 
  paste0(din_legen, "/dereplicated/orthogroups_occurrence.csv")

# create output  folder 
if (! dir.exists(dout)) dir.create(dout)

# read genome metadata table
genomes <- read_csv(fin_genomes, col_types = cols())

# read table with per-species occurrence estimations of orthogroups
orthogroup_occurrences <- 
  read_csv(fin_orthogroup_occurrences, col_types = cols())

###################
# EXPLORE GENOMES # 
###################

genomes %>%
  filter(quality >= 0.90) %>%
  count(species, name = "n_genomes") %>%
  ggplot(aes(x = rank(- n_genomes, ties.method = "first"), y = n_genomes)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 100, lty = 2) + 
  scale_y_log10() +
  theme_bw() +
  xlab("species") + ylab("high-quality genomes")
ggsave(
  paste0(dout, "/genome_counts.png"), units = "cm", width = 14, 
  height = 6
)

#####################
# EXPLORE PANGENOME # 
#####################

# determine species occurrence and core/accessory occurrence of orthogroups
orthogroups <-
  orthogroup_occurrences %>%
  mutate(status = if_else(occurrence_est > 0.8, "core", "accessory")) %>%
  group_by(orthogroup) %>%
  summarize(
    n_species = n(),
    n_species_core = sum(status == "core"),
    n_species_accessory = sum(status == "accessory")
  )

# visualize core/accessory counts of orthogroups in scatterplot 
n_species <- max(orthogroups$n_species)
orthogroups %>%
  ggplot(aes(x = n_species_core, y = n_species_accessory)) +
  geom_point() +
  xlim(c(0, n_species)) + ylim(c(0, n_species)) + 
  theme_bw() +
  xlab("number of species where core") + 
  ylab("number of species where accessory")
ggsave(
  paste0(dout, "/core_vs_accessory_counts.png"), units = "cm", width = 14, 
  height = 14
)

# visualize accessory status for orthogroups present in max 10 species
orthogroups %>%
  filter(n_species < 10, n_species > 1) %>%
  count(n_species, n_species_accessory, name = "n_orthogroups_accessory") %>%
  ggplot(aes(x = as.factor(n_species_accessory), y = n_orthogroups_accessory)) +
  geom_col() +
  facet_grid(~ n_species, scales = "free", space = "free_x") +
  theme_bw() +
  xlab("species where accessory") +
  ylab("number of orthogroups")
  # scale_y_log10() 
ggsave(
  paste0(dout, "/accessory_counts_max9_species.png"), units = "cm", width = 16, 
  height = 12
)

# visualize fixation frequency for orthogroups present in 10 or more species
orthogroups_filtered <- 
  orthogroups %>%
  filter(n_species >= 10) %>%
  mutate(fixation_frequency = n_species_core / n_species)
fig_acc_counts <- 
  orthogroups_filtered %>%
  ggplot(aes(x = fixation_frequency, fill = cut_number(n_species, 5))) +
  geom_histogram() +
  theme_bw() +
  xlab("fixation frequency (%)") +
  ylab("number of orthogroups") +
  scale_fill_brewer(palette = "OrRd", name = "number of species")
fig_acc_counts
ggsave(
  paste0(dout, "/fixation_frequency_min10_species.png"), units = "cm", width = 16, 
  height = 12
)
saveRDS(fig_acc_counts, paste0(dout, "/fixation_frequency_min10_species.rds"))

# inspect summary statistics of accessory status
orthogroups_filtered %>%
  summarize(
    n_species_min10genomes = length(unique(orthogroup_occurrences$species)), 
    n_orthogroups_min10species = n(), 
    p_orthogroups_max_5p_ff = sum(fixation_frequency < 0.05) / n(), 
    p_orthogroups_min_95p_ff = sum(fixation_frequency >= 0.95) / n()
  ) %>%
  write_csv(paste0(dout, "/accessory_analysis_statistics.csv"))

# write table that identifies most accessory orthogroups
orthogroups_filtered %>%
  filter(fixation_frequency < 0.05) %>%
  arrange(desc(n_species)) %>%
  select(orthogroup, n_species, n_species_accessory) %>%
  write_csv(paste0(dout, "/orthogroups_systematically_accessory.csv"))

# explore species-unique core orthogroups 
# (for Inas' species-species primers paper)
signatures <- 
  orthogroups %>%
  filter(n_species_core == 1, n_species_accessory == 0) %>%
  pull(orthogroup)
species_genomes <- 
  genomes %>%
  count(species, name = "n_genomes")
species <- 
  orthogroup_occurrences %>%
  filter(orthogroup %in% {{signatures}}) %>%
  count(species, name = "n_signatures") %>%
  right_join(species_genomes, by = "species") %>%
  filter(n_genomes >= 10) %>%
  replace_na(list(n_signatures = 0)) %>%
  arrange(n_signatures)
species %>%
  write_csv(paste0(dout, "/species_signature_genes.csv"))
