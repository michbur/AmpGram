options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 50
))

my_DT <- function(x, ...)
  datatable(x, ..., escape = FALSE, extensions = 'Buttons', filter = "top", rownames = FALSE,
            style = "bootstrap")

AMP_DT <- function(x, ...) {
  df <- x
  colnames(df) <- c("Putative AMP", "Probability")
  formatRound(my_DT(df, ...), 2, 4) 
}


plot_single_protein <- function(single_prot) {
  p <- ggplot(single_prot, aes(x = start, xend = end,
                               y = pred, yend = pred, color = decision,
                               linetype = decision)) +
    geom_segment() +
    geom_hline(yintercept = 0.5, color = "red") +
    ggtitle(single_prot[["seq_name"]][1]) +
    scale_x_continuous("Position") +
    scale_y_continuous("Probability of AMP", limits = c(0, 1)) +
    scale_color_manual("AMP", values = c(No = "#878787", Yes = "black")) + 
    scale_linetype_manual("AMP", values = c(No = "dashed", Yes = "solid")) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  if(max(single_prot[["end"]] > 100))
    p <- p + scale_x_continuous("Position", breaks = seq(0, max(single_prot[["end"]]), 
                                                         by = 20))
  
  p
}

predict_in_shiny <- function(object, newdata) {
  
  ngrams <- object[["imp_features"]]
  
  decoded_ngrams <- gsub(pattern = "_", replacement = ".", 
                         x = biogram::decode_ngrams(ngrams), fixed = TRUE)
  prediction_percentage <- 0
  withProgress(message = "", value = 0, {
    all_preds <- lapply(1L:length(newdata), function(ith_seq_id) {
      ith_seq <- toupper(newdata[[ith_seq_id]])
      ngram_count <- AmpGram:::find_ngrams(seq = ith_seq, decoded_ngrams = decoded_ngrams)
      colnames(ngram_count) <- ngrams
      all_mers_pred <- predict(object[["rf_mers"]], ngram_count)[["predictions"]][, 2]
      single_prot_pred <- predict(object[["rf_peptides"]], 
                                  AmpGram:::calculate_statistics(all_mers_pred))[["predictions"]][, 2]
      res <- list(seq = ith_seq,
                  all_mers_pred = all_mers_pred,
                  single_prot_pred = single_prot_pred)
      
      class(res) <- "single_ampgram_pred"
      
      prediction_percentage <<- prediction_percentage + 1/length(newdata)*100
      incProgress(1/length(newdata), detail = paste0(round(prediction_percentage, 2), 
                                                     "% proteins analyzed"))
      
      res
    })
    
  }, style = "old")
  
  if(is.null(names(newdata))) {
    names(all_preds) <- paste0("seq", 1L:length(all_preds))
  } else {
    names(all_preds) <- names(newdata)
  }
    
  all_preds
}
