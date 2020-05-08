library(shiny)
library(ggplot2)
library(AmpGram)
library(AmpGramModel)
library(DT)
library(shinythemes)
library(shinycssloaders)
library(markdown)

data(AmpGram_model)

options(shiny.maxRequestSize=10*1024^2)

options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 50
))

my_DT <- function(x, ...)
  datatable(x, ..., escape = FALSE, extensions = 'Buttons', filter = "top", rownames = FALSE,
            style = "bootstrap")

plot_single_protein <- function(single_prot) {
  p <- ggplot(single_prot, aes(x = start, xend = end,
                               y = pred, yend = pred, color = decision,
                               linetype = decision)) +
    geom_segment() +
    geom_hline(yintercept = 0.5, color = "red") +
    ggtitle(single_prot[["seq_name"]][1]) +
    scale_x_continuous("Position") +
    scale_y_continuous("Probability of AMP", limits = c(0, 1)) +
    scale_color_manual("AMP", values = c("#878787", "black")) + 
    scale_linetype_manual("AMP", values = c("dashed", "solid")) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  if(max(single_prot[["end"]] > 100))
    p <- p + scale_x_continuous("Position", breaks = seq(0, max(single_prot[["end"]]), 
                                                         by = 20))
  
  p
}

shinyServer(function(input, output) {
  
  prediction <- reactive({
    
    if (!is.null(input[["seq_file"]]))
      input_sequences <- read_txt(input[["seq_file"]][["datapath"]])
    input[["use_area"]]
    isolate({
      if (!is.null(input[["text_area"]]))
        if(input[["text_area"]] != "")
          input_sequences <- read_txt(textConnection(input[["text_area"]]))
    })
    
    if(exists("input_sequences")) {
      if(length(input_sequences) > 50) {
        #dummy error, just to stop further processing
        stop("Too many sequences. Please use AmpGram locally.")
      } else {
        if(any(lengths(input_sequences) < 10)) {
          #dummy error, just to stop further processing
          stop("The minimum length of the sequence is 10 amino acids.")
        } else {
          predict(AmpGram_model, input_sequences)
        }
      }
    } else {
      NULL
    }
  })
  
  decision_table <- reactive({
    if(!is.null(prediction())) {
      pred2df(prediction())
    }
  })
  
  output[["decision_table"]] <- renderDataTable({
    df <- decision_table()
    colnames(df) <- c("Protein name", "AMP probability")
    my_DT(df) %>% 
      formatRound(2, 4) 
    
  })
  
  detailed_preds <- reactive(({
    validate(
      need(input[["decision_table_rows_selected"]], 
           "Select at least one row in the Results table")
    )
    
    selected_pred_data <- prediction()[input[["decision_table_rows_selected"]]]
    
    detailed_pred_list <- lapply(1L:length(selected_pred_data), function(ith_pred_id) {
      ith_pred <- selected_pred_data[[ith_pred_id]]
      
      data.frame(seq_name = names(selected_pred_data)[ith_pred_id],
                 start = 1L:length(ith_pred[["all_mers_pred"]]), 
                 end = 1L:length(ith_pred[["all_mers_pred"]]) + 9, 
                 pred = ith_pred[["all_mers_pred"]],
                 decision = ith_pred[["all_mers_pred"]] > 0.5)
    })
    
  }))
  
  
  output[["detailed_preds"]] <- renderUI({
    detailed_preds_list <- lapply(1L:length(detailed_preds()), function(i) {
      list(plotOutput(paste0("detailed_plot", i)))
    })
    c(list(downloadButton("download_long", "Download long output (without graphics)"),
           downloadButton("download_long_graph", "Download long output (with graphics)")),
      do.call(tagList, unlist(detailed_preds_list, recursive = FALSE)))
  })
  
  
  for (i in 1L:300) {
    local({
      my_i <- i
      
      output[[paste0("detailed_plot", my_i)]] <- renderPlot(plot_single_protein(detailed_preds()[[my_i]]))
    })
  }
  
  output[["detailed_tab"]] <- renderUI({
    uiOutput("detailed_preds")
  })
  
  output[["dynamic_ui"]] <- renderUI({
    if (!is.null(input[["seq_file"]]))
      input_sequences <- read_txt(input[["seq_file"]][["datapath"]])
    input[["use_area"]]
    isolate({
      if (!is.null(input[["text_area"]]))
        if(input[["text_area"]] != "")
          input_sequences <- read_txt(textConnection(input[["text_area"]]))
    })
    
    if(exists("input_sequences")) {
      fluidRow(
        tags$p(HTML("<h3><A HREF=\"javascript:history.go(0)\">Start a new query</A></h3>"))
      )
    } else {
      fluidRow(
        h4("Exemplary sequences"),
        pre(includeText("prots.txt"))
      )
    }
  })
  
  output[["dynamic_tabset"]] <- renderUI({
    if(is.null(prediction())) {
      
      tabPanel(title = "Sequence input",
               tags$textarea(id = "text_area", style = "width:90%",
                             placeholder="Paste sequences (FASTA format required) here...", rows = 22, cols = 60, ""),
               p(""),
               actionButton("use_area", "Submit data from field above"),
               p(""),
               fileInput('seq_file', 'Submit .fasta or .txt file:'))
      
      
    } else {
      tabsetPanel(
        tabPanel("Results (tabular)",
                 dataTableOutput("decision_table")
        ),
        tabPanel("Results (graphical)",
                 uiOutput("detailed_tab")
        )
      )
    }
  })
  
  file_name <- reactive({
    if(is.null(input[["seq_file"]][["name"]])) {
      part_name <- "AmpGram_results"
    } else {
      part_name <- strsplit(input[["seq_file"]][["name"]], ".", fixed = TRUE)[[1]][1]
    }
    part_name
  })
  
  output[["download_long"]] <- downloadHandler(
    filename  = function() { 
      paste0(file_name(), "_pred.txt") 
    },
    content <- function(file) {
      sink(file, type = "output")
      cat("Input file name: ", ifelse(is.null(input[["seq_file"]][["name"]]), "none",
                                      input[["seq_file"]][["name"]]), "\n\n")
      cat(paste0("Date: ", Sys.time()), "\n\n")
      for (i in 1L:length(prediction())) {
        cat("\n\n")
        summary(prediction()[[i]])
        cat("\n\n")
      }
      sink()
    }
  )
  
  output[["download_long_graph"]] <- downloadHandler(
    filename  = function() { 
      paste0(file_name(), "_pred.html") 
    },
    content <- function(file) {
      src <- normalizePath("AmpGram-report.Rmd")
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, "AmpGram-report.Rmd")
      out <- render("AmpGram-report.Rmd", output_format = "html_document", file, quiet = TRUE)
      file.rename(out, file)
    }
  )
  
})
