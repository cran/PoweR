help.law <- function(law.index) {

    # Retrieve the number of laws in our package:
    tmp <- names(getDLLRegisteredRoutines("PoweR")[[".C"]])
    nb.laws <- length(grep("law", tmp))
  
    if (!(law.index %in% 1:nb.laws)) stop(paste("Law index should be an integer between 1 and ", nb.laws, ".", sep = ""))
  
    if (nchar(law.index) == 1) Rd <- paste0("law000", law.index)
    if (nchar(law.index) == 2) Rd <- paste0("law00", law.index)
    if (nchar(law.index) == 3) Rd <- paste0("law0", law.index)
    
    help(Rd)
    
}
