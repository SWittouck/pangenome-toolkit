#' Calculate pangenome statistics for a single pangenome
#' 
#' @param pan A tibble with the columns gene, genome and orthogroup. 
#' @param orthogroup A quoted alternative name for the orthogroup column. 
#' 
#' @return A named list with pangenome statistics
pangenome_stats_one <- function(pan, orthogroup = "orthogroup") {
  
  statistics <- list()
  
  # preparation
  pan <- pan %>% select(gene, genome, orthogroup = all_of(orthogroup))
  orthogroups <- 
    pan %>%
    count(genome, orthogroup, name = "n_copies") %>%
    group_by(orthogroup) %>%
    summarize(
      n_genomes = n(),
      n_genomes_sc = sum(n_copies == 1),
      n_genes = sum(n_copies),
      average_cn = n_genes / n_genomes,
      .groups = "drop"
    )
  
  # number of orthogroups and genomes
  statistics$n_orthogroups <- nrow(orthogroups)
  n_genomes <- pan$genome %>% unique() %>% length()
  
  # number of core orthogroups
  statistics$n_core <- sum(orthogroups$n_genomes == n_genomes)
  
  # number of single-copy core orthogroups
  statistics$n_hard_scg <- sum(orthogroups$n_genomes_sc == n_genomes)
  
  # number of soft single-copy core orthogroups
  statistics$n_soft_scg <- sum(orthogroups$n_genomes_sc >= 0.95 * n_genomes)
  
  # genes in single-copy orthogroups
  orthogroups_sc <- orthogroups[orthogroups$average_cn == 1, ][["orthogroup"]]
  statistics$perc_in_single_copy <-
    pan %>%
    filter(orthogroup %in% {{orthogroups_sc}}) %>%
    {nrow(.) / nrow(pan)}
  
  statistics
  
}

#' Calculate pangenome statistics for all pangenomes in a table
#' 
#' @param pans A tibble with the columns gene and genome, and one column
#'   per pangenome tool with the orthogroups inferred by that tool. 
#'   
#' @return A tibble with the column tool and a bunch of columns with pangenome
#'   statistics. 
pangenome_stats_all <- function(pans) {
  
  names(pans) %>%
    setdiff(c("gene", "genome")) %>%
    {structure(., names = .)} %>%
    map(~ pangenome_stats_one(pans, orthogroup = .x)) %>%
    map2(names(.), ~ {.x$tool <- .y; .x}) %>%
    transpose() %>%
    map(as_vector) %>%
    as_tibble() %>%
    relocate(tool)
  
}

#' Compile a table with benchmark measures for various tools on a single dataset
#' 
#' @param din A folder with the files "pangenomes.csv" and "resources.csv"
#' 
#' @return A tibble with the column "tool" and a bunch of other columns
compile_tooltable <- function(din) {
  
  # define read_csv variant that reads all cols as character
  read_csv_c <- function(...) read_csv(..., col_types = cols(.default = "c"))
  
  # read pangenome and resource consumption tables
  pans <- read_csv(paste0(din, "/pangenomes.csv"), col_types = cols()) 
  resources <- read_csv_c(paste0(din, "/resources.csv")) 
  
  # compute pangenome statistics
  tools_panstats <- pans %>% pangenome_stats_all() 
  
  # join all tool benchmarks and return
  tools_panstats %>% full_join(resources, by = "tool")
  
}

benchmark_plot <- 
  function(tools, toolnames, toollabels, ncol = 3, measures = "all") {
  
  # convert hms time to minutes
  tools$elapsed_time_mins <-
    tools$elapsed_time %>%
    if_else(str_count(., ":") == 1, str_c("0:", .), .) %>%
    lubridate::hms() %>% 
    lubridate::period_to_seconds() %>%
    {. / 60}
  
  # define the order and labels of the statistics to plot
  stats <- c(
    "n_orthogroups", "n_core", "elapsed_time_mins", "f_score", "precision", 
    "recall", "n_hard_scg", "n_soft_scg", "perc_in_single_copy"
  )
  statlabels <- c(
    "# families (x1000)", "# core families", "time (mins)", "F-score", 
    "precision", "recall", "# SCOs", "# soft SCOs", "% in single-copy family"
  )
  
  # prepare tooltable for figure 
  tools <-
    tools %>%
    select(tool, any_of(stats)) %>%
    mutate(n_orthogroups = n_orthogroups / 1000) %>%
    pivot_longer(cols = - tool, names_to = "measure", values_to = "value") %>%
    mutate(na_label = if_else(is.na(value), "NA", as.character(NA))) %>%
    mutate(tool = factor(tool, levels = toolnames, labels = toollabels)) %>%
    mutate(measure = factor(measure, levels = stats, labels = statlabels)) %>%
    mutate(tool = fct_rev(tool))
  
  if (measures[1] != "all") tools <- tools %>% filter(measure %in% measures)
  
  # make the actual figure
  tools %>%
    ggplot(aes(x = tool, y = value)) +
    geom_col() + 
    geom_text(aes(x = tool, y = 0, label = na_label), hjust = 0) +
    facet_wrap(~ measure, scales = "free_x", ncol = ncol) +
    xlab("") + ylab("") +
    coord_flip() +
    theme_classic()
  
}

#' Compute the f-measure
#' 
#' Given a vector with predicted groups and a vector with reference groups of 
#' objects, this function calculates the precision, recall and f-measure of the
#' predicted groups.  
#' 
#' @param group_pred A vector with predicted groups. 
#' @param group_ref A vector with reference groups. 
#' 
#' @return A list with the elements "precision", "recall" and "f_measure".
f_measure <- function(group_pred, group_ref) {
  
  pan <- tibble(cluster_pred = group_pred, cluster_ref = group_ref)
  
  counts <- 
    pan %>%
    count(cluster_ref, cluster_pred)
  clusters_ref <-
    counts %>%
    group_by(cluster_ref) %>%
    arrange(desc(n)) %>%
    slice(1) %>%
    ungroup() %>%
    rename(majority_pred_cluster = cluster_pred) %>%
    select(- n)  
  clusters_pred <-
    counts %>%
    group_by(cluster_pred) %>%
    arrange(desc(n)) %>%
    slice(1) %>%
    ungroup() %>%
    rename(majority_ref_cluster = cluster_ref) %>%
    select(- n)
  pan %>%
    left_join(clusters_ref, by = "cluster_ref") %>%
    left_join(clusters_pred, by = "cluster_pred") %>%
    summarize(
      precision = sum(majority_ref_cluster == cluster_ref) / n(),
      recall = sum(majority_pred_cluster == cluster_pred) / n()
    ) %>%
    {list(
      precision = .$precision, recall = .$recall,
      f_measure = 2 / ((1 / .$precision) + (1 / .$recall))
    )}
  
}
