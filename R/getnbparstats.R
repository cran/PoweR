getnbparstats <- function(stat.indices = NULL) {
    
    
    tmp <- names(getDLLRegisteredRoutines("PoweR")[[".C"]])
    ind.stats <- as.numeric(sub("stat", "", tmp[grep("stat[1-9]", tmp)]))
    if (!all(stat.indices %in% ind.stats)) stop(paste("The values", paste(stat.indices[which(!(stat.indices %in% ind.stats))], collapse = " "),
                                                      " in 'stat.indices' do not correspond to any defined statistic!", collapse = ""))
    
    if (is.null(stat.indices)) stat.indices <- ind.stats
    
    lst <- length(stat.indices)
    nbparstats.list <- as.vector(rep(0, lst))
    
  for (i in 1:lst) {
        Cstat.name <- "tmp" # To remove a NOTE at R CMD check
        Cstat.name <- paste("stat", stat.indices[i], sep = "")
        #statsym <- getNativeSymbolInfo(Cstat.name, PACKAGE = "PoweR")
mydotC <- get(".PoweR_stat_dispatch", envir = asNamespace("PoweR"))[[Cstat.name]]; if (is.null(mydotC)) stop("Unknown stat function: ", Cstat.name)

        nbparstats.list[i] <- (#.C(statsym, 
        mydotC(0.0, 1L, 0.0, 1L, rep(" ", 50), 1L, 0.0, 0L,
                                  0.0, 0.0, 0.0, 0L, 0L, 0L, 0.0, nbparamstat = 0L, PACKAGE = "PoweR"))$nbparamstat
    }
    
    return(nbparstats.list)
    
}
