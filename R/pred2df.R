#' @param x results of prediction as produced by \code{\link{predict.ampgram_model}}
#' @export
pred2df <- function(x) {
  data.frame(seq_name = names(x),
             probability = sapply(x, function(i) i[["single_prot_pred"]]))
}
