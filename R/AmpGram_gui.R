#' AmpGram Graphical User Interface
#'
#' Launches graphical user interface that predicts presence of 
#' antimicrobial peptides.
#'
#' @importFrom shiny runApp
#' @section Warning : Any ad-blocking software may cause malfunctions.
#' @export AmpGram_gui
AmpGram_gui <- function() {
  require_AmpGramModel()
  runApp(system.file("AmpGram", package = "AmpGram"))
}