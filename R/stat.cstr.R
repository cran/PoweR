stat.cstr <- function(stat.index, stat.pars = NULL, n = 0) {


  # We get the default (maximal) number of parameters
  Cstat.name <- "tmp" # To remove a NOTE at R CMD check
  Cstat.name <- paste("stat", stat.index, sep = "")
  # We get the total number of parameters
  #statsym <- getNativeSymbolInfo(Cstat.name, PACKAGE = "PoweR")
mydotC <- get(".PoweR_stat_dispatch", envir = asNamespace("PoweR"))[[Cstat.name]]; if (is.null(mydotC)) stop("Unknown stat function: ", Cstat.name)

  out <- #.C(statsym, 
  mydotC(0.0, as.integer(n), 0.0, 0L, rep(" ", 50), 1L, 0.0, 0L, 0.0, 0.0, 0.0, 0L, 0L, 0L, 0.0, nbparams = 1L, PACKAGE = "PoweR")
  nbparams <- out$nbparams
  if (length(stat.pars) > nbparams) stop(paste("Length of 'stat.pars' should be at most", nbparams))
  # We get the default values of the parameters (using the trick of putting the first value of *name to "1"
  out <- #.C(statsym, 
  mydotC(0.0, as.integer(n), 0.0, 0L, name = c("1", rep(" ", 49)), 1L, 0.0, 0L, 0.0, 0.0, 0.0,
            0L, alter = 0L, 0L, params = c(stat.pars, rep(0.0, nbparams - length(stat.pars))), nbparams = as.integer(nbparams), NAOK = TRUE, PACKAGE = "PoweR")
  if (length(stat.pars) >= 1) {
    if (n == 0) {
      stat.pars <- c(stat.pars, out$params[-(1:length(stat.pars))])
    } else { # the default value of some parameters might depend on n
      stat.pars <- out$params[1:nbparams]
    }
  }
  if (is.null(stat.pars)) stat.pars <- out$params[1:nbparams]
    
  name <- out$name
  name <- gsub('\\', '', gsub('$', '', sub(' +$', '', paste(out$name, collapse = "")), fixed = TRUE), fixed = TRUE)
  
  return(list(name = name, nbparams = nbparams, stat.pars = stat.pars, alter = if (out$alter == 0) c("0,1,2") else out$alter)) 
}
