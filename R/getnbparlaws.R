getnbparlaws <- function(law.indices = NULL) {
  
  
    tmp <- names(getDLLRegisteredRoutines("PoweR")[[".C"]])
    ind.laws <- grep("law", tmp[grep("law",tmp)])
    if (!all(law.indices %in% ind.laws)) stop(paste("The values", paste(law.indices[which(!(law.indices %in% ind.laws))], collapse = " "),
                                                    " in 'law.indices' do not correspond to any defined law!", collapse = ""))
    
    if (is.null(law.indices)) law.indices <- ind.laws
    
    lst <- length(law.indices)
    nbparlaws.list <- as.vector(rep(0, lst))
    
    for (i in 1:lst) {
      Claw.name <- paste("law", law.indices[i], sep = "")
      #lawsym <- getNativeSymbolInfo(Claw.name, PACKAGE = "PoweR")
mydotC <- get(".PoweR_law_dispatch", envir = asNamespace("PoweR"))[[Claw.name]]; if (is.null(mydotC)) stop("Unknown law function: ", Claw.name)

        nbparlaws.list[i] <- (#.C(lawsym, 
        mydotC(1L, 0.0, rep(" ", 50), 1L, rep(0.0, 4),
                                 nbparamlaw = 1L, 1L, PACKAGE = "PoweR"))$nbparamlaw
    }
    
    return(nbparlaws.list)
    
}
