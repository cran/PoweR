create.alter <- function(stat.indices = c(42,51,61), values.alter = NULL) {

  if(getRversion() < "3.1") dontCheck <- identity

  if (is.null(stat.indices) || any(is.na(stat.indices))) {
    stop("stat.indices should not contain NULL or NA values.")
  }
  
  nbstats <- length(stat.indices)
  
  alter <- as.list(c())
  
  for (i in 1:nbstats) {
    
    stat.index <- stat.indices[i]
    value.alter <- values.alter[i]

    if (is.null(value.alter) || is.na(value.alter)) {
      # We retrieve the default value:
      alter[[paste("stat",stat.index,sep="")]] <- (.C(dontCheck(paste("stat",as.character(stat.index),sep="")),0.0,1L,0.0,1L,rep(" ",50),1L,0.0,0L,0.0,0.0,0.0,0L,0L,0L,0.0,0L,PACKAGE="PoweR"))[[13]]
    } else {
      alter[[paste("stat",stat.index,sep="")]] <- value.alter
    }
  }
  
  return(alter)
}
