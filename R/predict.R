#' Predict antimicrobial peptides
#'
#' Recognizes antimicrobial peptides using the AmpGram algorithm.
#' @param object \code{ampgram_model} object.
#' @param newdata \code{list} of sequences (for example as given by
#' \code{\link[biogram]{read_fasta}}).
#' @param ... further arguments passed to or from other methods.
#' @export
#' @importFrom biogram binarize decode_ngrams
#' @importFrom pbapply pblapply
#' @importFrom ranger ranger
#' @importFrom stats predict
#' @importFrom stringi stri_count

predict.ampgram_model <- function(object, newdata, ...) {
  require_AmpGramModel()
  
  ngrams <- object[["imp_features"]]
    
  decoded_ngrams <- gsub(pattern = "_", replacement = ".", 
                         x = decode_ngrams(ngrams), fixed = TRUE)
  
  
  all_preds <- pblapply(newdata, function(ith_seq) {
    ngram_count <- find_ngrams(seq = ith_seq, decoded_ngrams = decoded_ngrams)
    colnames(ngram_count) <- ngrams
    all_mers_pred <- predict(object[["rf_mers"]], ngram_count)[["predictions"]][, 2]
    single_prot_pred <- predict(object[["rf_peptides"]], 
                                calculate_statistics(all_mers_pred))[["predictions"]][, 2]
    res <- list(seq = ith_seq,
                all_mers_pred = all_mers_pred,
                single_prot_pred = single_prot_pred)
    
    class(res) <- "single_ampgram_pred"
    
    res
  })
  
  if(is.null(names(all_preds))) 
    names(all_preds) <- paste0("seq", 1L:length(all_preds))
  
  all_preds
}


# data(AmpGram_model)
# sample_seq <- list(seq1 = c("F", "E", "N", "C", "N", "I", "T", "M", "G", "N", "M", "V", 
#                             "R", "H", "I", "R", "W", "Y", "R", "D", "R", "Q", "K", "G", "D", 
#                             "Y", "W", "W", "Y", "T", "I", "K", "Y", "S", "M", "A", "M", "I", 
#                             "A", "C", "N", "I", "N", "V", "T", "I", "N", "Q", "C", "V"),
#                    seq2 = c("Q", "Y", "T", "S", "I", "M", "F", "L", "T", "A", "G", "H", 
#                             "L", "A", "P", "W", "D", "R", "W", "C", "R", "S", "L", "T", "T", 
#                             "W", "F", "G", "A", "P", "S", "A", "T", "Y", "P", "F", "F", "W", 
#                             "E", "P", "E", "D", "I", "I", "I", "K", "P", "N", "T", "A"))
# predict(AmpGram_model, sample_seq)
