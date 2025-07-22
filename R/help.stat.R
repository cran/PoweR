help.stat <- function(stat.index) {

    # Retrieve the number of stats in our package:
    tmp <- names(getDLLRegisteredRoutines("PoweR")[[".C"]])
    nb.stats <- length(grep("stat", tmp))

    if (!(stat.index %in% 1:nb.stats)) stop(paste("Statistic index should be an integer between 1 and ", nb.stats, ".", sep = ""))
  
  
    if (nchar(stat.index) == 1) Rd <- paste0("stat000", stat.index)
    if (nchar(stat.index) == 2) Rd <- paste0("stat00", stat.index)
    if (nchar(stat.index) == 3) Rd <- paste0("stat0", stat.index)
    
    help(Rd)
    
}
