many.crit <- function(law.index,stat.indices,M=10^3,vectn=c(20,50,100),levels=c(0.05,0.1),alter=create.alter(stat.indices),law.pars=NULL,parstats=NULL,model=NULL) {

  if(getRversion() < "3.1") dontCheck <- identity
  
  stats.len <- length(stat.indices)

# Management of alter
  if (!is.null(alter)) {
    if (!is.list(alter)) stop("'alter' should be a list")
    if (is.null(names(alter))) stop("'parstats' should be a named list")
    if (any(is.na(names(alter)))) stop("'alter' names should all be defined")
    if (length(alter) != length(stat.indices)) stop("'alter' and 'stat.indices' should have the same length")
    for (s in 1:stats.len) {
      if (names(alter)[s] != paste("stat",stat.indices[s],sep="")) stop(paste("Name of 'alter'[[",s,"]] should be equal to 'stat",stat.indices[s],sep=""))
      if (!(alter[[s]] %in% 0:4)) stop(paste("'alter'[[",s,"]] should be in  {0,1,2,3,4}.",sep=""))
      Cstat.name <- paste("stat",as.character(stat.indices[s]),sep="")
      alter.true <- .C(dontCheck(Cstat.name),as.double(0.0),1L,0.05,1L,rep(" ",50),1L,0.0,0L,0.0,0.0,0.0,0L,alter=as.integer(alter[s]),0L,rep(0.0,4),0L,PACKAGE="PoweR")$alter
      if (alter[[s]] != alter.true) {
        warning(paste("'alter'[[",s,"]] should be set to ",alter.true,". We have done this for you!"),sep="")
        alter[[s]] <- alter.true
      }
    }
    alter <- unlist(alter)
  } else { # alter is NULL
    alter <- rep(0,stats.len)
    names(alter) <- paste("stat",stat.indices,sep="")
  }

 
# Management of parstats
 nbparstats <- getnbparstats(stat.indices)
 if (!is.null(parstats)) {
   if (!is.list(parstats)) stop("'parstats' should be a list")
   if (is.null(names(parstats))) stop("'parstats' should be a named list")
   if (any(is.na(names(parstats)))) stop("'parstats' names should all be defined")
   if (length(parstats) != length(stat.indices)) stop("'parstats' and 'stat.indices' should have the same length")
   for (s in 1:stats.len) {
     if (names(parstats)[s] != paste("stat",stat.indices[s],sep="")) stop(paste("Name of 'parstats'[[",s,"]] should be equal to 'stat",stat.indices[s],sep=""))
     if (!is.na(parstats[[s]]) && (nbparstats[s] == 0)) stop(paste("'parstats[['",s,"]] should be equal to NA",sep=""))
     if ((nbparstats[s] != 0) && (length(parstats[[s]]) != nbparstats[s])) stop(paste("The length of parstats[[",s,"]] should be ",nbparstats[s],sep=""))
   }
 }

 
 nblevels <- length(levels)
 vectn.len <- length(vectn)
 mylist <- as.list(c())
 for (s in 1:stats.len) {
 
   stat.index <- stat.indices[s]
   statname <- paste("stat",as.character(stat.index),sep="")
    
   mylist[[s]] <- matrix(NaN,nrow=length(vectn)*length(levels),ncol=2,dimnames=list(NULL,c("critL","critR")))
   
   for (l in 1:nblevels) {  
     for (n in 1:vectn.len) {
       if (alter[[s]] == 0) { # two.sided test
         mylist[[s]][n+vectn.len*(l-1),] <- compquant(n=vectn[n],law.index=law.index,stat.index=stat.index,probs=c(levels[l]/2,1-levels[l]/2),M=M,law.pars=law.pars,stat.pars=parstats[[s]],model=model)$quant
       } else if ((alter[[s]] == 1) || (alter[[s]] == 4)) { # less test, or bilateral test that reject H0 only for small values of the test statistic
         res <- compquant(n=vectn[n],law.index=law.index,stat.index=stat.index,probs=levels[l],M=M,law.pars=law.pars,stat.pars=parstats[[s]],model=model)$quant
         mylist[[s]][n+vectn.len*(l-1),] <- c(res,NA)
       } else if ((alter[[s]] == 2) || (alter[[s]] == 3)) { # greater test, or bilateral test that reject H0 only for large values of the test statistic
         res <- compquant(n=vectn[n],law.index=law.index,stat.index=stat.index,probs=1-levels[l],M=M,law.pars=law.pars,stat.pars=parstats[[s]],model=model)$quant
         mylist[[s]][n+vectn.len*(l-1),] <- c(NA,res)
       }
     }
   }

   # add a new column named "param" to distinguish the critical values from the same statistic with different parameters
   # if parstats is NULL, "param" will appear as NA
   mat <- expand.grid(vectn,levels,paste(parstats[[s]],collapse=" "))
   colnames(mat) <- c("n","level","param")
   mylist[[s]] <- cbind(mat,mylist[[s]])

   names(mylist)[s] <- statname
   
 }

 list.names.replicated <- names(mylist)[table(names(mylist)) > 1]
 
 for (name in list.names.replicated) {

   mask <- names(mylist) %in% name

   names(mylist)[mask] <- paste(name,1:sum(mask),sep=".")
      
 }
  
 class(mylist) <- "critvalues"
 
 return(mylist)
 
 
}





























## First solution for creating vector t :
# nb <- c(2,2,3)
# t <- 1:sum(nb)
# toto <- c()
# for (s in 1:length(nb)) {
  
  # toto <- c(toto,rep(s,nb[s]))

# }
# t <- split(t,factor(toto))


## Second solution for creating vector t : 
# apply(cbind(cumsum(x) - x + 1,x),MARGIN=1,FUN=function(x){seq(from=x[1],length=x[2])})
   
# t <- apply(cbind(cumsum(nbparstats2) - nbparstats2 + 1,nbparstats2),MARGIN=1,FUN=function(x){seq(from=x[1],length=x[2])})
 
 
 
