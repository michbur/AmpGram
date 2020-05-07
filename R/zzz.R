.onLoad <- function(libname, pkgname) {
  options(AmpGram_suppress_prompt = FALSE)
  if ("AmpGramModel" %in% rownames(installed.packages()))
    warning("To be able to use AmpGram properly, you should install 'AmpGramModel' package available via GitHub. You can do it by calling 'install_AmpGramModel()'")
}
