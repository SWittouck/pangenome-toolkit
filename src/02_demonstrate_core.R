# This script evaluates the SCARAP core module, by making a comparison with the
# full pangenome on a dataset with L. plantarum genomes.
# 
# Abbreviations: 
# - coremod = core module 
# - panmod = pan module

# dependencies: R v4.2.3, tidyverse v1.3.1, eulerr v6.1.1

library(tidyverse)
source("src/02_demonstrate_core_functions.R")

#################################################
# Read data and process into one big gene table #
#################################################

# define relevant paths 
fin_core <- "data/benchmark_scarap_v3/scarap_core/lactiplantibacillus/corefull/genes.tsv"
fin_pan <- "data/benchmark_scarap_v3/scarap_core/lactiplantibacillus/pan/pangenome.tsv"
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

###################################
# Create one big orthogroup table #
###################################

# calculate single-copy occurrence of pan module orthogroups
orthogroups_panmod <-
  genes %>%
  group_by(orthogroup_panmod) %>%
  summarize(
    genes_panmod = n(), occurrence_panmod = sum(table(genome) == 1)
  )

# calculate single-copy occurrence of core module orthogroups
orthogroups_coremod <-
  genes %>%
  filter(! is.na(orthogroup_coremod)) %>%
  group_by(orthogroup_coremod) %>%
  summarize(
    genes_coremod = n(), occurrence_coremod = sum(table(genome) == 1)
  )

# map coremod orthogroups to panmod orthogroups with greatest gene content
# overlap
orthogroups_coremod2 <-
  genes %>%
  filter(! is.na(orthogroup_coremod)) %>%
  count(orthogroup_coremod, orthogroup_panmod, name = "shared_genes") %>%
  arrange(desc(shared_genes)) %>%
  group_by(orthogroup_coremod) %>%
  slice(1) %>%
  ungroup()

# make one big orthogroup table
orthogroups <- 
  orthogroups_panmod %>% 
  left_join(orthogroups_coremod2, by = "orthogroup_panmod") %>%
  left_join(orthogroups_coremod, by = "orthogroup_coremod") %>%
  select(
    orthogroup_panmod, orthogroup_coremod, occurrence_panmod, 
    occurrence_coremod, genes_coremod, genes_panmod, shared_genes
  )

# which panmod orthogroups have more than one coremod orthogroup mapped to them? 
orthogroups %>%
  filter(occurrence_panmod >= 0.95 * max(occurrence_panmod)) %>%
  group_by(orthogroup_panmod) %>%
  summarize(n_core_orthogroups = sum(! is.na(orthogroup_coremod))) %>%
  filter(n_core_orthogroups != 1)

# retain only one coremod orthogroup per panmod orthogroup
orthogroups <- 
  orthogroups %>%
  group_by(orthogroup_panmod) %>%
  arrange(desc(shared_genes)) %>%
  slice(1) %>%
  ungroup() 

# remove intermediate datasets
rm(orthogroups_panmod, orthogroups_coremod, orthogroups_coremod2)

#########################################################
# Compare the core genome from the core and pan modules #
#########################################################

# calculate the number of genomes
n_genomes <- length(genes$genome %>% unique())

# average percentage shared genes of core orthogroups? 
orthogroups %>%
  filter(! is.na(orthogroup_coremod)) %>%
  {mean(.$shared_genes / .$genes_coremod)}

# what is the gene-level precision and recall of the core module?
coregenes_panmod <- 
  orthogroups %>%
  filter(occurrence_panmod >= 0.95 * n_genomes) %>%
  pull(orthogroup_panmod)
statistics_genelevel <- 
  genes %>%
  mutate(core_panmod = orthogroup_panmod %in% coregenes_panmod) %>%
  mutate(core_coremod = ! is.na(orthogroup_coremod)) %>%
  {precrec(.$core_coremod, .$core_panmod)} %>%
  mutate(level = "genes")

# what is the orthogroup-level precision and recall of the core module? 
statistics_orthogrouplevel <- 
  orthogroups %>%
  mutate(core_panmod = occurrence_panmod >= 0.95 * n_genomes) %>%
  mutate(core_coremod = ! is.na(orthogroup_coremod)) %>%
  {precrec(.$core_coremod, .$core_panmod)} %>%
  mutate(level = "orthogroups")

# write statistics 
bind_rows(statistics_genelevel, statistics_orthogrouplevel) %>%
  relocate(level) %>%
  write_csv(paste0(dout, "/statistics.csv"))

# REMARK: the rest of the script is currently not used in the manuscript!! 

#############################################################
# Visualize the comparison between the core and pan modules #
#############################################################

# show the gene-level overlap in a venn diagram
fout <- paste0(dout, "/genes_venn.png")
png(fout, units = "cm", width = 8, height = 6, res = 400)
statistics_genelevel %>%
  {c(`coremod&panmod` = .$tp, `coremod` = .$fp, `panmod` = .$fn)} %>%
  eulerr::euler(shape = "ellipse") %>%
  plot(
    quantities = T, lty = 0, fills = c("#a6cee3", "#1f78b4"), 
    labels = c("core\nmodule", "pan\nmodule")
  )
dev.off() 

# show the orthogroup-level overlap in a venn diagram
fout <- paste0(dout, "/orthogroups_venn.png")
png(fout, units = "cm", width = 8, height = 6, res = 400)
statistics_orthogrouplevel %>%
  {c(`coremod&panmod` = .$tp, `coremod` = .$fp, `panmod` = .$fn)} %>%
  eulerr::euler(shape = "ellipse") %>%
  plot(
    quantities = T, lty = 0, fills = c("#a6cee3", "#1f78b4"), 
    labels = c("core\nmodule", "pan\nmodule")
  )
dev.off() 

# visualize the comparison between the core and pangenome for a random subsample
# of 100 orthogroups (with high occurrence)
orthogrouplvls <-
  orthogroups %>%
  arrange(desc(occurrence_panmod)) %>%
  pull(orthogroup_panmod)
orthogroups_sub <- 
  orthogroups %>%
  filter(occurrence_panmod > 100) %>%
  sample_n(100) %>%
  mutate(occurrence_panmod = occurrence_panmod / {{n_genomes}}) %>%
  mutate(occurrence_coremod = occurrence_coremod / {{n_genomes}})
orthogroups_sub %>%
  mutate(
    orthogroup_panmod = factor(orthogroup_panmod, levels = orthogrouplvls)
  ) %>%
  ggplot(aes(x = orthogroup_panmod, y = occurrence_panmod)) +
  geom_col() +
  geom_point(aes(y = occurrence_coremod)) +
  geom_hline(yintercept = 0.95) + 
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("pan module orthogroup") + ylab("occurrence")
ggsave(
  paste0(dout, "/coremod_vs_panmod_randomsample.png"), units = "cm", 
  width = 16, height = 10
)
