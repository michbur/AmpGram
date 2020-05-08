#' Convert predictions to data.frame
#' Return predictions as data.frame
#' @param x results of prediction as produced by \code{\link{predict.ampgram_model}}
#' @return a data.frame containing sequences with their probabilities 
#' @export
pred2df <- function(x) {
  data.frame(seq_name = names(x),
             probability = sapply(x, function(i) i[["single_prot_pred"]]))
}
