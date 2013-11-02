powcomp.fast <- function(law.indices,stat.indices,vectn=c(20,50,100),M=10^3,levels=c(0.05,0.1),critval=NULL,alter=create.alter(stat.indices),parlaws=NULL,parstats=NULL,nbclus=1,model=NULL,null.law.index=2,null.law.pars=NULL) {
  
  if(getRversion() < "3.1") dontCheck <- identity

  if (nbclus>1) {
#    suppressWarnings(parallel.pkg.present <- require(parallel))
    parallel.pkg.present <- is.loaded("parallel")
    suppressWarnings(Rmpi.pkg.present <- is.loaded("Rmpi"))
    if (all(!c(parallel.pkg.present,Rmpi.pkg.present))) stop("Either package parallel or Rmpi should be installed!")
#    suppressWarnings(rsprng.pkg.present <- require(rsprng))
#    if (!rsprng.pkg.present) stop("Package rsprng is not installed!")
    cluster.type <- if (parallel.pkg.present) "PSOCK" else "MPI" # We prefer to use "PSOCK" (i.e. parallel) because it's easier.
  }
  
  vectn.len <- length(vectn)
  stats.len <- length(stat.indices)
  laws.len  <- length(law.indices)
  nblevel   <- length(levels)  
  
# Management of critval and creation of critvalL, critvalR  and usecrit
# If we provide a single value in critval$statj then it is critvalR
# If we provide two values in critval$statj then it is c(critvalL,critvalR)  ... IN THAT ORDER!!
  critvalL <- critvalR <- rep(0,vectn.len*stats.len*nblevel)
  usecrit <- rep(0,vectn.len*stats.len*nblevel)
  if (is.null(critval)) {
    warning(paste("'critval' has been computed internally using function many.crit() with the value of 'law.index'=",null.law.index," (i.e. ",law.cstr(null.law.index)$name,")",sep=""))
    critval <- many.crit(law.index=null.law.index,stat.indices,M,vectn,levels,alter,null.law.pars,parstats)
  }
  if (!is.list(critval)) stop("'critval' should be a list")
  if (is.null(names(critval))) stop("'critval' should be a named list")
  if (any(is.na(names(critval)))) stop("'critval' names should all be defined")
  if (length(critval) != length(stat.indices)) stop("'critval' and 'stat.indices' should have the same length")
  for (s in 1:stats.len) {
    crittmp <- as.matrix(critval[[s]][,-3])
    if (ncol(crittmp) == 1) crittmp <- t(crittmp)
    for (l in 1:nblevel) {
      vals <- crittmp[crittmp[,"level"]==levels[l],]
      for (k in 1:vectn.len) {
        if (vectn.len == 1) {
          vals2 <- vals[-(1:2)]
        } else {
          vals2 <- vals[vals[,"n"]==vectn[k],-(1:2)]
        }
        if (!is.null(vals2)) {
          usecrit[k+vectn.len*(l-1)+nblevel*vectn.len*(s-1)] <- 1
          critvalL[k+vectn.len*(l-1)+nblevel*vectn.len*(s-1)] <- vals2[1]
          critvalR[k+vectn.len*(l-1)+nblevel*vectn.len*(s-1)] <- vals2[2]
        }
      }      
    }
  }
  
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
 

# Management of parlaws
  parlawtmp <- nbparlaws <- c()
  if (!is.null(parlaws)) {
    if (!is.list(parlaws)) stop("'parlaws' should be a list")
    if (is.null(names(parlaws))) stop("'parlaws' should be a named list")
    if (any(is.na(names(parlaws)))) stop("'parlaws' names should all be defined")
    if (length(parlaws) != length(law.indices)) stop("'parlaws' and 'law.indices' should have the same length")
    for (s in 1:laws.len) {
      if (names(parlaws)[s] != paste("law",law.indices[s],sep="")) stop(paste("Name of 'parlaws'[[",s,"]] should be equal to 'law",law.indices[s],sep=""))
      if (length(parlaws[[s]]) > 4) stop(paste("Length of 'parlaws'[[",s,"]] should not exceed 4.",sep=""))
      if ((length(parlaws[[s]]) > 1) && any(is.na(parlaws[[s]]))) stop(paste("'parlaws'[[",s,"]] cannot contain NA values unless its length is 1.",sep=""))
      nbparlaws <- c(nbparlaws,length(na.omit(parlaws[[s]])))
      parlawtmp <- c(parlawtmp,c(parlaws[[s]],rep(0.0,4-nbparlaws[s])))
    }
  } else {
    for (s in 1:laws.len) {
      tmp <- law.cstr(law.indices[s])
      nbparlaws <- c(nbparlaws,tmp$nbparams)
      parlawtmp <- c(parlawtmp,c(tmp$law.pars,rep(0.0,4-nbparlaws[s])))
    }
  }
  parlaws <- parlawtmp

  
# Management of parstats
  nbparstats <- getnbparstats(stat.indices)
  parstatstmp <- c()
  if (!is.null(parstats)) {
    if (!is.list(parstats)) stop("'parstats' should be a list")
    if (is.null(names(parstats))) stop("'parstats' should be a named list")
    if (any(is.na(names(parstats)))) stop("'parstats' names should all be defined")
    if (length(parstats) != length(stat.indices)) stop("'parstats' and 'stat.indices' should have the same length")
    for (s in 1:stats.len) {
      if (names(parstats)[s] != paste("stat",stat.indices[s],sep="")) stop(paste("Name of 'parstats'[[",s,"]] should be equal to 'stat",stat.indices[s],sep=""))
      if (!is.na(parstats[[s]]) && (nbparstats[s] == 0)) stop(paste("'parstats[['",s,"]] should be equal to NA",sep=""))
      if ((nbparstats[s] != 0) && (length(parstats[[s]]) != nbparstats[s])) stop(paste("The length of parstats[[",s,"]] should be ",nbparstats[s],sep=""))
      parstatstmp <- c(parstatstmp,parstats[[s]])
    }
  } else {
    for (s in 1:stats.len) parstatstmp <- c(parstatstmp,stat.cstr(stat.indices[s])$stat.pars)
  }
  parstats <- parstatstmp
  parstats[is.na(parstats)] <- 0


# Management of model
  ## FAUDRA REGARDER CA DE PLUS PRES QUAND JE RENDRAI MODEL FONCTIONNEL!!  
  if (is.double(model) || is.integer(model)) {
    modelnum <- model
    funclist <- list(function(){})
    thetavec <- 0
    xvec <- 0
    p <- length(thetavec)
    np <- length(xvec)
  } else {
    if (is.null(model)) {
      modelnum <- 1
      funclist <- list(function(){})
      thetavec <- 0
      xvec <- 0
      p <- length(thetavec)
      np <- length(xvec)
    } else { # model should be a list (function(x,thetavec,xvec),theta,xvec)
      modelnum <- 0
      funclist <- list(model[[1]])
      thetavec <- model[[2]]
      xvec <- model[[3]]
      p <- length(thetavec)
      np <- length(xvec)     
    }
  }



# We perform the computations

  decision.len <- stats.len*vectn.len*laws.len*nblevel
  decision <- rep(0,decision.len)

  # a) Using a cluster
  if (nbclus > 1) { # We start the cluster
    # makeCluster = Create a set of copies of R running in parallel and communicating over sockets or using MPI.
    cl <- makeCluster(nbclus, type = cluster.type)		
#    clusterSetupSPRNG(cl)
                                        
      
    myfunc <- function(M) {
      require(PoweR)
      .C("powcompfast",M=as.integer(M),law.indices=as.integer(law.indices),laws.len=as.integer(laws.len),vectn=as.integer(vectn),vectn.len=as.integer(vectn.len),stat.indices=as.integer(stat.indices),stats.len=as.integer(stats.len),decision=as.integer(decision),decision.len=as.integer(decision.len),levels=as.double(levels),nblevel=as.integer(nblevel),
         cL=as.double(critvalL),cR=as.double(critvalR),usecrit=as.integer(usecrit),alter=as.integer(alter),nbparlaws=as.integer(nbparlaws),parlaws=as.double(parlaws),nbparstats=as.integer(nbparstats),parstats=as.double(parstats),as.integer(modelnum), funclist, as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),PACKAGE="PoweR",NAOK=TRUE)
    }
    
    out <- clusterCall(cl, myfunc, round(M/nbclus)) # M/nbclus iterations are performed on each core
          
      # We stop the cluster
    stopCluster(cl)
    
  } else {

  #b) or without a cluster

    out <- list(.C("powcompfast",M=as.integer(M),law.indices=as.integer(law.indices),laws.len=as.integer(laws.len),vectn=as.integer(vectn),vectn.len=as.integer(vectn.len),stat.indices=as.integer(stat.indices),stats.len=as.integer(stats.len),decision=as.integer(decision),decision.len=as.integer(decision.len),levels=as.double(levels),nblevel=as.integer(nblevel),
                  cL=as.double(critvalL),cR=as.double(critvalR),usecrit=as.integer(usecrit),alter=as.integer(alter),nbparlaws=as.integer(nbparlaws),parlaws=as.double(parlaws),nbparstats=as.integer(nbparstats),parstats=as.double(parstats),as.integer(modelnum), funclist, as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),PACKAGE="PoweR",NAOK=TRUE))

  }


    out[[1]] <- out[[1]][c("M","law.indices","vectn","stat.indices","decision","levels","cL","cR","usecrit","alter","nbparlaws","parlaws","nbparstats","parstats")]

  if (nbclus > 1) {
    for (clus in 2:nbclus) {
      out[[1]]$decision <- out[[1]]$decision + out[[clus]]$decision
    }
  }

  out[[1]]$nbclus <- nbclus
  
  k <- 1
  for (i in 1:length(out[[1]]$nbparstats)) {
    if (out[[1]]$nbparstats[i] == 0) {out[[1]]$parstats[k]<- NA ; k <- k + 1} else k <- k + out[[1]]$nbparstats[i]
  }

  return(structure(out[[1]], class = c("power","list")))
  
}

