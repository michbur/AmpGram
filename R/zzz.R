.onLoad <- function(libname, pkgname) {
  options(AmpGram_suppress_prompt = FALSE)
  check_AmpGramModel <- try(find.package("AmpGramModel"), silent = TRUE)
  if (inherits(check_AmpGramModel, "try-error"))
    warning("To be able to use AmpGram properly, you should install 'AmpGramModel' package available via GitHub. You can do it by calling 'install_AmpGramModel()'")
}
