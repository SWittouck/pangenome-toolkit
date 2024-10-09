#' Compile a table with benchmark measures for various tools on a single dataset
#' 
#' @param din A folder with the files "statistics.csv", "benchmarks.csv" 
#'   (optional) and "resources.csv". 
#' 
#' @return A tibble
compile_runtable <- function(din) {
  
  # define read_csv variant that reads all cols as character
  read_csv_c <- function(...) read_csv(..., col_types = cols(.default = "c"))

  # define paths 
  fin_stats <- paste0(din, "/statistics.csv")
  fin_bench <- paste0(din, "/benchmarks.csv")
  fin_resources <- paste0(din, "/resources.csv")

  # read tables
  runs_stats <- read_csv(fin_stats, col_types = cols())
  if (file.exists(fin_bench)) 
    runs_bench <- read_csv(fin_bench, col_types = cols()) 
  runs_resources <- read_csv_c(fin_resources)

  # merge tables
  runs <- runs_stats 
  if (file.exists(fin_bench)) runs <- full_join(runs, runs_bench, by = "tool")
  runs <- full_join(runs, runs_resources, by = "tool")

  # return
  runs %>% rename(run = tool)
  
}

#' Create a benchmark barplot
#' 
#' @param runtable A table with measurements of pangenome runs; should contain 
#'   the columns "run", "dataset" and var. 
#' @param var Unquoted variable to plot. 
#' @param xlab Label to plot on the x-axis. 
#' 
#' @return A facet wrapped bar plot (ggplot)
runplot_bar <- function(runtable, var, xlab) {
  var <- enquo(var)
  runtable %>%
    ggplot(aes(x = run, y = {{var}})) +
    geom_col() + 
    geom_label(aes(y = 0, label = {{var}}), hjust = 0, size = 3) +
    facet_wrap(facets = "dataset", scales = "free_x") +
    xlab("") + ylab(xlab) +
    coord_flip() +
    theme_bw()
}

#' Create a benchmark time plot
#' 
#' @param runtable A table with measurements of pangenome runs; should contain 
#'   the columns "run", "dataset" and "elapsed_time". 
#' @param ytext Logical value indicated whether run names should be shown on the
#'   y-axis or not. 
#' 
#' @return A facet wrapped dot plot (ggplot)
runplot_time <- function(runtable, ytext = T) {
  plot <- 
    runtable %>%
    mutate(
      time_mins = elapsed_time %>%
        if_else(str_count(., ":") == 1, str_c("0:", .), .) %>%
        lubridate::hms() %>% 
        lubridate::period_to_seconds() %>%
        {. / 60}
    ) %>%
    ggplot(aes(x = run, y = time_mins)) +
    geom_point(shape = "|", size = 3) + 
    facet_wrap(facets = "dataset") +
    scale_y_log10() + 
    xlab("") + ylab("time (minutes)") +
    coord_flip() +
    theme_bw() 
  if (! ytext) {
    no_y <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    plot <- plot + no_y
  }
  plot
}

#' Add letter to plot
#' 
#' This function adds a capital letter to a plot for panel arranging purposes. 
#' 
#' @param plot A ggplot2 object
#' @param letter A character value
#' 
#' @return An updated ggplot2 object
addletter <- function(plot, letter) {
  g <- ggplotGrob(plot + ggtitle(letter))
  g$layout$l[g$layout$name == "title"] <- 1
  g
}
