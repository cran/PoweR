statcompute <- function(stat.index, data, levels = c(0.05, 0.1),
                        critvalL = NULL, critvalR = NULL, alter = 0,
                        stat.pars = NULL, pvalcomp = 1L, check = TRUE) {

  n <- length(data)
  nblevels <- length(levels)
  Cstat.name <- "tmp" # To remove a NOTE at R CMD check
  Cstat.name <- paste("stat", stat.index, sep = "")

  cL <- 0    
  cR <- 0    

  usecrit <- 1

  if (check) {  # The following instructions take time! This is why the 'check'
                # argument exists.
    tmp <- names(getDLLRegisteredRoutines("PoweR")[[".C"]])
    ind.stats <- grep("stat", tmp)
    nb.stats <- length(ind.stats)
    if (!(paste("stat", stat.index, sep = "") %in% tmp))
      stop("This test statistic has not been included in the package!")

    nbparstat <- getnbparstats(stat.index)
    if (is.null(stat.pars) || is.na(stat.pars)) {
      stat.pars <- rep(0.0, nbparstat) # C++ technical requirement.
      nbparstat <- 0 # The default values will be used by the C++ function.
    } else {
      tmp <- length(stat.pars)
      stat.pars <- c(stat.pars, rep(0.0, nbparstat - length(stat.pars)))
      nbparstat <- tmp
    }
    if (n < 2) stop("'data' should be a vector of length at least 2")
    
    if (!(alter %in% 0:4)) stop("'alter' should be 0, 1, 2, 3 or 4")
#    statsym <- getNativeSymbolInfo(Cstat.name, PACKAGE = "PoweR")
mydotC <- get(".PoweR_stat_dispatch", envir = asNamespace("PoweR"))[[Cstat.name]]; if (is.null(mydotC)) stop("Unknown stat function: ", Cstat.name)

    testname <- #.C(statsym, 
    mydotC(0.0, 0L, 0.0, 0L,
                   statname = rep(" ",50), 1L, 0.0, 0L, 0.0, 0.0, 0.0, 0L,
                   alter = 0L, 0L, rep(0.0, 4), 1L, PACKAGE = "PoweR")$statname
    testname <- gsub('\\', '',
                     gsub('$', '',
                          sub(' +$', '', paste(testname, collapse = "")),
                          fixed = TRUE), fixed = TRUE) # Remove $ and backlash signs
    if (!is.null(critvalL)) {
      if (length(critvalL) != nblevels)
        stop("'critvalL' length and 'levels' length should not differ!")
      cL <- critvalL
    }
    
    if (!is.null(critvalR)) {
      if (length(critvalR) != nblevels)
        stop("'critvalR' length and 'levels' length should not differ!")
      cR <- critvalR
    }

    if (is.null(critvalL) && is.null(critvalR)) usecrit <- 0

  } else {
    nbparstat <- length(stat.pars)
    testname <- ''
  }
#  statsym <- getNativeSymbolInfo(Cstat.name, PACKAGE = "PoweR")
mydotC <- get(".PoweR_stat_dispatch", envir = asNamespace("PoweR"))[[Cstat.name]]; if (is.null(mydotC)) stop("Unknown stat function: ", Cstat.name)

  out <- #.C(statsym, 
  mydotC(as.double(data), as.integer(n),
            as.double(levels), as.integer(nblevels), rep(" ", 50),
            0L, statistic = 0.0, pvalcomp = pvalcomp, pvalue = 0.0,
            cL = as.double(cL), cR = as.double(cR), as.integer(usecrit),
            alter = as.integer(alter), decision = as.integer(rep(0, nblevels)),
            stat.pars = as.double(stat.pars), nbparstat = as.integer(nbparstat),
            PACKAGE = "PoweR")
  
  if (out$pvalcomp == 0L) out$pvalue <- NA
  
  return(list(statistic = out$statistic, pvalue = out$pvalue,
              decision = out$decision, alter = out$alter,
              stat.pars = out$stat.pars[1:out$nbparstat], symbol = testname))
  
}

