#' @export
#' @importFrom devtools install_github
install_AmpGramModel <- function() {
  install_github("michbur/AmpGramModel")
}

require_AmpGramModel <- function() {
  if (!suppressWarnings(require("AmpGramModel", quietly = TRUE)) && !getOption("AmpGram_suppress_prompt")) {
    response <- menu(c("yes", "no", "no and don't ask me anymore"), 
                     title = "To be able to use AmpGram properly, you should have installed 'AmpGramModel' package available via GitHub. Install?")
    switch (response, 
            tryCatch(install_AmpGramModel(),
               finally = if (!suppressWarnings(require("AmpGramModel"))) 
                 warning("There was an error during an attempt to install 'AmpGramModel' package.", call. = FALSE)),
            warning("You cannot access full functionality of this package without having installed 'AmpGramModel'. You can do it manually by calling 'devtools::install_github('michbur/AmpGramModel')'", call. = FALSE),
            {options(AmpGram_suppress_prompt = TRUE)
              cat("Ok, but you cannot access full functionality of this package without having installed 'AmpGramModel'")},
            warning("You cannot access full functionality of this package without having installed 'AmpGramModel'. You can do it manually by calling 'devtools::install_github('michbur/AmpGramModel')'", call. = FALSE)
            )
  }
}


