#!/usr/bin/env Rscript

# This script explores the benchmarking results of the SCARAP pan module and
# other pangenome tools on various datasets.

# dependencies: R v4.2.3, tidyverse v1.3.1, ggpubr v0.4.0

library(tidyverse)
source("src/01_demonstrate_pan_functions.R")

# define input/output folders
din <- "data/benchmark_scarap_v4/scarap_pan"
dout <- "results/pan_demonstration"

# create output  folder 
if (! dir.exists(dout)) dir.create(dout)

#############################
# Prepare benchmarking data #
#############################

# read data from scarap pan runs 
runtables <- 
  list.files(din) %>%
  paste0(din, "/", .) %>%
  set_names(str_extract(., "[^/]+$")) %>%
  map(compile_runtable)

# prepare proper names and order of datasets and runs/tools
runlabels <-
  c(
    "orthofinder_blast" = "OrthoFinder-BLAST", 
    "orthofinder_mmseqs" = "OrthoFinder-MMseqs2", 
    "sonicparanoid_sensitive" = "SonicParanoid-sensitive", 
    "sonicparanoid_fast" = "SonicParanoid-fast", 
    "broccoli" = "Broccoli",  
    "pirate" = "PIRATE",
    "scarap_s" = "SCARAP-no-splits",
    "scarap_fh" = "SCARAP"
  )
datalabels <- 
  c(
    "repseqs" = "repseqs",
    "simpan" = "SimPan",
    "tonkinhill" = "Tonkin-Hill",
    "lactobacillales_genusreps" = "lacto-genusreps",
    "orthobench" = "OrthoBench",
    "parabench" = "paraBench"
  )

# extract and preprocess the repseq runs
repseqruns <- 
  runtables$repseqs %>%
  mutate(f_measure = round(f_measure, 3)) %>%
  mutate(dataset = "lacto-species")
runtables$repeqs < NULL

# preprocess the non-repseq runtables
runtables$lactobacillales_genusreps$f_measure <- NA
runtables <- 
  runtables %>%
  map2(names(.), function(rt, dataset) {
    rt %>% 
      mutate(run = {{runlabels}}[run]) %>%
      mutate(run = factor(run, levels = {{runlabels}})) %>%
      mutate(run = fct_rev(run)) %>%
      mutate(dataset = {{datalabels}}[{{dataset}}]) %>%
      mutate(dataset = factor(dataset, levels = {{datalabels}})) %>%
      mutate(f_measure = round(f_measure, 3))
  })

##########################
# Explore the benchmarks #
##########################

# visualize min_reps benchmarks
repseqs_min_reps <- 
  repseqruns %>%
  filter(str_detect(run, "maxreps0p_maxalign1f$")) %>%
  mutate(run = str_extract(run, "^[^_]+")) %>%
  mutate(run = fct_reorder(run, as.integer(str_extract(run, "[0-9]+")))) %>%
  mutate(run = fct_rev(run))
ggpubr::ggarrange(
  repseqs_min_reps %>%
    runplot_bar(f_measure, "f-measure"),
  repseqs_min_reps %>%
    runplot_time(ytext = F),
  ncol = 2, nrow = 1, widths = c(3, 2)
)
ggsave(
  paste0(dout, "/min_reps.png"), units = "cm", width = 16, height = 6
)

# visualize max_align benchmarks
repseqs_max_align <- 
  repseqruns %>%
  filter(str_detect(run, "^minreps32_maxreps0p")) %>%
  mutate(run = str_extract(run, "[^_]+$")) %>%
  mutate(run = fct_reorder(run, as.integer(str_extract(run, "[0-9]+")))) %>%
  mutate(run = fct_rev(run))
ggpubr::ggarrange(
  repseqs_max_align %>%
    runplot_bar(f_measure, "f-measure"),
  repseqs_max_align %>%
    runplot_time(ytext = F),
  ncol = 2, nrow = 1, widths = c(3, 2)
)
ggsave(
  paste0(dout, "/max_align.png"), units = "cm", width = 16, height = 6
)

# visualize benchmarks of prokaryotic datasets 
runtables2 <- map(runtables, ~ filter(., run != "SCARAP-no-splits"))
ggpubr::ggarrange(
  runtables2$simpan %>% 
    runplot_bar(f_measure, "f-measure") %>% 
    addletter("A"),
  runtables2$simpan %>% 
    runplot_time(ytext = F) %>% 
    addletter("B"),
  runtables2$tonkinhill %>% 
    runplot_bar(f_measure, "f-measure") %>% 
    addletter("C"),
  runtables2$tonkinhill %>% 
    runplot_time(ytext = F) %>% 
    addletter("D"),
  runtables2$lacto %>% 
    runplot_bar(sc_core_orthogroups, "# SCOs") %>% 
    addletter("E"),
  runtables2$lacto %>% 
    runplot_time(ytext = F) %>% 
    addletter("F"),
  ncol = 2, nrow = 3, widths = c(3, 2)
)
ggsave(
  paste0(dout, "/prok_datasets.png"), units = "cm", width = 16, height = 20
)

# visualize benchmarks on eukaryotic datasets
ggpubr::ggarrange(
  runtables2$orthobench %>% 
    runplot_bar(f_measure, "f-measure") %>% 
    addletter("A"),
  runtables2$orthobench %>% 
    runplot_time(ytext = F) %>% 
    addletter("B"),
  runtables2$parabench %>% 
    runplot_bar(f_measure, "f-measure") %>% 
    addletter("C"),
  runtables2$parabench %>% 
    runplot_time(ytext = F) %>% 
    addletter("D"),
  ncol = 2, nrow = 2, widths = c(3, 2)
)
ggsave(
  paste0(dout, "/euk_datasets.png"), units = "cm", width = 16, height = 14
)

# compare scarap-s to scarap-nf 
scarap_runs <- 
  runtables %>%
  reduce(bind_rows) %>% 
  filter(dataset != "repseqs", str_detect(run, "^SCARAP")) 
ggpubr::ggarrange(
  scarap_runs %>% 
    runplot_bar(sc_core_orthogroups, "single-copy core orthogroups") %>% 
    addletter("A"),
  scarap_runs %>% 
    filter(dataset != "lacto-genusreps") %>%
    runplot_bar(f_measure, "f-measure") %>% 
    addletter("B"),
  ncol = 1, nrow = 2
)
ggsave(
  paste0(dout, "/scarap-no-split.png"), units = "cm", width = 16, height = 14
)
