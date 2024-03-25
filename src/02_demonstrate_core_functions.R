#' Compute the precision and recall
#' 
#' Given a boolean vector with predictions and a boolean vector with the ground
#' truth, return the number of true positives, false positives, false negatives
#' and true negaties, as well as the precision and recall. 
#' 
#' @param group_pred A vector with predicted groups. 
#' @param group_ref A vector with reference groups. 
#' 
#' @return A list with the elements "precision", "recall" and "f_measure".
precrec <- function(pred, truth) {
  tibble(pred, truth) %>%
    summarize(
      tp = sum(pred & truth),
      fn = sum(! pred & truth),
      fp = sum(pred & ! truth),
      tn = sum(! pred & ! truth),
      precision = tp / (tp + fp),
      recall = tp / (tp + fn)
    ) 
}