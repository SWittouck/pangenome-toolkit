# This script explores the benchmarking results of the SCARAP pan module and
# other pangenome tools on various datasets.

# dependencies: R v4.2.3, tidyverse v1.3.1, ggpubr v0.4.0

library(tidyverse)
source("src/01_demonstrate_pan_functions.R")

#######################################################################
# Read the benchmarking data and format into a nice table per dataset #
#######################################################################

# define input/output folders
din <- "data/benchmark_scarap_v3/scarap_pan"
dout <- "results/pan_demonstration"

# define paths of input subfiles/subfolders
din_genusreps <- paste0(din, "/lactobacillales_genusreps")
din_speciesreps <- paste0(din, "/lactobacillales_speciesreps")
din_orthobench <- paste0(din, "/orthobench")
din_parabench <- paste0(din, "/parabench")

# create output  folder 
if (! dir.exists(dout)) dir.create(dout)

# read and process pangenome tables and add resource consumption tables
tools_genusreps <- compile_tooltable(din_genusreps)
tools_orthobench <- compile_tooltable(din_orthobench)
tools_parabench <- compile_tooltable(din_parabench)

# read orthobench and parabench benchmark tables 
benchmarks_orthobench <- 
  paste0(din_orthobench, "/benchmarks.csv") %>%
  read_csv(col_types = cols()) %>%
  mutate(across(c(f_score, precision, recall), ~ parse_number(.) / 100))
benchmarks_parabench <- 
  paste0(din_parabench, "/benchmarks.csv") %>%
  read_csv(col_types = cols()) %>%
  rename(f_score = f1_score)

# add orthobench benchmarks to orthobench tooltable (if not yet present)
if (! "f_score" %in% names(tools_orthobench)) {
  tools_orthobench <- 
    tools_orthobench %>%
    full_join(benchmarks_orthobench, by = "tool")
}

# add parabench benchmarks to parabench tooltable (if not yet present)
if (! "f_score" %in% names(tools_parabench)) {
  tools_parabench <- 
    tools_parabench %>%
    full_join(benchmarks_parabench, by = "tool")
}

# define tool names (for internal use) and tool labels (for figures/tables)
toolnames <-
  c(
    "orthofinder_blast", "orthofinder_mmseqs", "sonicparanoid_sensitive", 
    "sonicparanoid_fast", "broccoli", "scarap_s", "scarap_fh"
  )
toollabels <- 
  c(
    "OrthoFinder-BLAST", "OrthoFinder-MMseqs2", "SonicParanoid-sensitive", 
    "SonicParanoid-fast", "Broccoli", "SCARAP-no-splits", "SCARAP"
  )

# set time of orthofinder_blast to zero because it's an outlier
tools_orthobench <-
  tools_orthobench %>%
  {.$elapsed_time[.$tool == "orthofinder_blast"] <- "NA"; .}
tools_parabench <-
  tools_parabench %>%
  {.$elapsed_time[.$tool == "orthofinder_blast"] <- "NA"; .}

# retain only the per-tool tables
rm(benchmarks_orthobench, benchmarks_parabench)

#########################################################
# Explore the benchmarks through figures and statistics #
#########################################################

# visualize benchmarks - summary plots
dout_sum <- paste0(dout, "/benchmarks_summary")
if (! dir.exists(dout_sum)) dir.create(dout_sum)
measures <- c("# families (x1000)", "# SCOs", "time (mins)")
(
  fig_genusreps <-
    tools_genusreps %>% 
    benchmark_plot(toolnames, toollabels, measures = measures)
)
ggsave(
  paste0(dout_sum, "/lactobacillales_genusreps.png"), units = "cm", width = 16, 
  height = 5
)
measures <- c("# families (x1000)", "F-score", "time (mins)")
(
  fig_orthobench <-
    tools_orthobench %>% 
    benchmark_plot(toolnames, toollabels, measures = measures)
)
ggsave(
  paste0(dout_sum, "/orthobench.png"), units = "cm", width = 16, height = 5
)
(
  fig_parabench <-
    tools_parabench %>% 
    benchmark_plot(toolnames, toollabels, measures = measures)
)
ggsave(
  paste0(dout_sum, "/parabench.png"), units = "cm", width = 16, height = 5
)

# make panel with summary plots 
give_letter <- function(plot, letter) {
  g <- ggplotGrob(plot + ggtitle(letter))
  g$layout$l[g$layout$name == "title"] <- 1
  g
}
ggpubr::ggarrange(
  give_letter(fig_genusreps, "A"),
  give_letter(fig_orthobench, "B"),
  give_letter(fig_parabench, "C"),
  ncol = 1
)
ggsave(
  paste0(dout_sum, "/panel.png"), units = "cm", width = 16, height = 18
)

# visualize benchmarks - extended plots
dout_ext <- paste0(dout, "/benchmarks_extended")
if (! dir.exists(dout_ext)) dir.create(dout_ext)
tools_genusreps %>% benchmark_plot(toolnames, toollabels)
ggsave(
  paste0(dout_ext, "/lactobacillales_genusreps.png"), units = "cm", width = 16, 
  height = 10
)
tools_orthobench %>% benchmark_plot(toolnames, toollabels)
ggsave(
  paste0(dout_ext, "/orthobench.png"), units = "cm", width = 16, height = 15
)
tools_parabench %>% benchmark_plot(toolnames, toollabels)
ggsave(
  paste0(dout_ext, "/parabench.png"), units = "cm", width = 16, height = 15
)

# visualize benchmarks - plots for slides
dout_slide <- paste0(dout, "/benchmarks_slide")
if (! dir.exists(dout_slide)) dir.create(dout_slide)
measures <- c("# families (x1000)", "# SCOs", "time (mins)")
tools_genusreps %>% benchmark_plot(toolnames, toollabels, measures = measures)
ggsave(
  paste0(dout_slide, "/lactobacillales_genusreps.png"), units = "cm", width = 16, 
  height = 6
)
measures <- c("F-score", "precision", "recall")
tools_orthobench %>% benchmark_plot(toolnames, toollabels, measures = measures)
ggsave(
  paste0(dout_slide, "/orthobench.png"), units = "cm", width = 16, height = 6
)
tools_parabench %>% benchmark_plot(toolnames, toollabels, measures = measures)
ggsave(
  paste0(dout_slide, "/parabench.png"), units = "cm", width = 16, height = 6
)

#############################################################################
# Compare the pangenome of the lacto genusreps with that of the speciesreps #
#############################################################################

# compare lacto_genusreps to lacto_speciesreps 
pan_genusreps <-
  paste0(din_genusreps, "/pangenomes.csv") %>%
  read_csv(col_types = cols()) %>%
  select(gene, genome, orthogroup_genusreps = scarap_fh)
colnames <- c("gene", "genome", "orthogroup_speciesreps")
pan_speciesreps <- 
  paste0(din_speciesreps, "/scarap_fh/pangenome.tsv") %>%
  read_tsv(col_names = colnames, col_types = cols())
pan_lacto <- 
  pan_genusreps %>%
  select(- genome) %>%
  left_join(pan_speciesreps, by = "gene")
c(
  f_measure(pan_lacto$orthogroup_speciesreps, pan_lacto$orthogroup_genusreps),
  pangenome_stats_one(pan_lacto, orthogroup = "orthogroup_speciesreps")
) %>%
  as_tibble() %>%
  pivot_longer(everything(), names_to = "measure", values_to = "value") %>%
  write_csv(paste0(dout, "/lactobacillales_speciesreps_subsetted.csv"))

# explore the differences between the pangenome of lacto_genusreps and
# lacto_speciesreps 
orthogroups_genusreps <- 
  pan_lacto %>%
  count(orthogroup_genusreps, orthogroup_speciesreps) %>%
  group_by(orthogroup_genusreps) %>%
  summarize(
    n_genes = sum(n),
    mode = max(n)
  )
orthogroups_genusreps %>%
  filter(n_genes > 1) %>%
  mutate(rank = rank(n_genes, ties.method = "first")) %>%
  pivot_longer(
    cols = c("n_genes", "mode"), names_to = "measure", values_to = "count"
  ) %>%
  ggplot(aes(x = rank, y = count, col = measure, group = measure)) +
  geom_point(size = 0.1) +
  scale_color_brewer(palette = "Paired") +
  theme_bw() 
ggsave(
  paste0(dout, "/genusreps_orthogroups_modes.png"), units = "cm", width = 16, 
  height = 12
)
orthogroups_genusreps %>%
  filter(n_genes > 1) %>%
  ggplot(aes(x = n_genes, y = mode)) +
  geom_jitter(size = 1) +
  # scale_x_log10() + 
  # scale_y_log10() + 
  theme_bw() 
ggsave(
  paste0(dout, "/genusreps_orthogroups_modes2.png"), units = "cm", width = 16, 
  height = 12
)
length(unique(pan_lacto$genome))
