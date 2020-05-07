.onLoad <- function(libname, pkgname) {
  options(AmpGram_suppress_prompt = FALSE)
  if (!suppressWarnings(require("AmpGramModel"))) 
    warning("To be able to use AmpGram properly, you should install 'AmpGramModel' package available via GitHub. You can do it by calling 'install_AmpGramModel()'")
}