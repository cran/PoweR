create.alter <- function(stat.indices = c(42, 51, 61), values.alter = NULL) {

    if (is.null(stat.indices) || any(is.na(stat.indices))) {
        stop("stat.indices should not contain NULL or NA values.")
    }
  
    nbstats <- length(stat.indices)
  
    alter <- as.list(c())
  
    for (i in 1:nbstats) {
    
        stat.index <- stat.indices[i]
        value.alter <- values.alter[i]
        
        Cstat.name <- "tmp" # To remove a NOTE at R CMD check
        Cstat.name <- paste("stat", as.character(stat.index), sep = "")

        if (is.null(value.alter) || is.na(value.alter)) {
          # We retrieve the default value:
          #statsym <- getNativeSymbolInfo(stat.name, PACKAGE = "PoweR")
          mydotC <- get(".PoweR_stat_dispatch", envir = asNamespace("PoweR"))[[Cstat.name]]; if (is.null(mydotC)) stop("Unknown stat function: ", Cstat.name)

            alter[[i]] <- (#.C(statsym, 
            mydotC(0.0, 1L, 0.0, 1L, rep(" ", 50), 1L, 0.0,
                                      0L, 0.0, 0.0, 0.0, 0L, 0L, 0L, 0.0, 0L, PACKAGE = "PoweR"))[[13]]
        } else {
            alter[[i]] <- value.alter
        }
      names(alter)[i] <- Cstat.name
    }
    
    return(alter)
}
