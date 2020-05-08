#' Get putative AMPs
#' @inheritParams pred2df
#' @export

get_AMPs <- function(x) {
  tenmer_start <- 1L:length(x[["all_mers_pred"]])
  only_AMP_start <- tenmer_start[x[["all_mers_pred"]] > 0.5]
  only_AMP_end <- only_AMP_start + 9
  data.frame(putative_AMP = sapply(1L:length(only_AMP_start), function(ith_pos) 
    paste0(x[["seq"]][only_AMP_start[ith_pos]:only_AMP_end[ith_pos]], collapse = "")
  ), prob = x[["all_mers_pred"]][x[["all_mers_pred"]] > 0.5])
}
