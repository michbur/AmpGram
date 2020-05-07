#' Predict amyloids
#'
#' Recognizes amyloids using AmpGram algorithm.
#' @param object \code{ampgram_model} object.
#' @param newdata \code{list} of sequences (for example as given by
#' \code{\link[biogram]{read_fasta}}).
#' @param ... further arguments passed to or from other methods.
#' @export
#' @importFrom biogram binarize decode_ngrams
#' @importFrom stringi stri_count
#' @examples
#' data(AmpGram_model)
#' sample_seq <- list(seq1 = c("F", "E", "N", "C", "N", "I", "T", "M", "G", "N", "M", "V", 
#'                             "R", "H", "I", "R", "W", "Y", "R", "D", "R", "Q", "K", "G", "D", 
#'                             "Y", "W", "W", "Y", "T", "I", "K", "Y", "S", "M", "A", "M", "I", 
#'                             "A", "C", "N", "I", "N", "V", "T", "I", "N", "Q", "C", "V"),
#'                    seq2 = c("Q", "Y", "T", "S", "I", "M", "F", "L", "T", "A", "G", "H", 
#'                             "L", "A", "P", "W", "D", "R", "W", "C", "R", "S", "L", "T", "T", 
#'                             "W", "F", "G", "A", "P", "S", "A", "T", "Y", "P", "F", "F", "W", 
#'                             "E", "P", "E", "D", "I", "I", "I", "K", "P", "N", "T", "A"))
#' predict(AmpGram_model, sample_seq)

predict.ampgram_model <- function(object, newdata, ...) {
  require_AmpGramModel()
  
  ngrams <- object[["imp_features"]]
    
  decoded_ngrams <- gsub(pattern = "_", replacement = ".", 
                         x = decode_ngrams(ngrams), fixed = TRUE)
  
  
  lapply(newdata, function(ith_seq) {
    ngram_count <- find_ngrams(seq = ith_seq, decoded_ngrams = decoded_ngrams)
    colnames(ngram_count) <- ngrams
    all_mers_pred <- predict(object[["rf_mers"]], ngram_count)[["predictions"]][, 2]
    single_prot_pred <- predict(object[["rf_peptides"]], 
                                calculate_statistics(all_mers_pred))[["predictions"]][, 2]
  })
  
}

count_longest <- function(x) {
  splitted_x <- strsplit(x = paste0(as.numeric(x > 0.5), collapse = ""),
                         split = "0")[[1]]
  len <- unname(sapply(splitted_x, nchar))
  if (length(len[len > 0]) == 0) {
    0 } else {
      len[len > 0]
    }
}

calculate_statistics <- function(pred) {
  data.frame(fraction_true = mean(pred > 0.5),
             pred_mean = mean(pred),
             pred_median = median(pred),
             n_peptide = length(pred),
             n_pos = sum(pred > 0.5),
             pred_min = min(pred),
             pred_max = max(pred), 
             longest_pos = max(count_longest(pred)),
             n_pos_10 = sum(count_longest(pred) >= 10),
             frac_0_0.2 = sum(pred <= 0.2)/length(pred),
             frac_0.2_0.4 = sum(pred > 0.2 & pred <= 0.4)/length(pred),
             frac_0.4_0.6 = sum(pred > 0.4 & pred <= 0.6)/length(pred),
             frac_0.6_0.8 = sum(pred > 0.6 & pred <= 0.8)/length(pred),
             frac_0.8_1 = sum(pred > 0.8 & pred <= 1)/length(pred)) 
}

find_ngrams <- function(seq, decoded_ngrams) {
  end_pos <- 10L:length(seq)
  start_pos <- end_pos - 9
  
  res <- binarize(do.call(rbind, lapply(1L:length(end_pos), function(ith_mer_id) {
    ten_mer <- paste0(seq[start_pos[ith_mer_id]:end_pos[ith_mer_id]], collapse = "")
    stri_count(ten_mer, regex = decoded_ngrams)
  })))
  
  res
}


