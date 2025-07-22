# Function in pure R to compute the criritcal values
many.critR <- function(law.H0, stat.indices, M = 10 ^ 3, vectn = c(20, 50, 100), levels = 0.05, alter,
                       law.pars = NULL, parstats) {

  eval(parse(text = paste("vecstats <- c(", paste("stat", stat.indices, sep = "", collapse = ", "), ")", sep = "")))

  for (i in 1:length(parstats)) {
    if (any(is.na(parstats[[i]]))) parstats[[i]] <- stat.cstr(as.numeric(strsplit(names(parstats[i]), "stat")[[1]][2]))$stat.pars
  }
  
  nmax <- max(vectn)
  nbstats <- length(vecstats)
  nbn <- length(vectn)
  nblevels <- length(levels)

  if (!is.null(law.pars)) {
    eval(parse(text = paste("law.H0 <- function(x, ", paste(law.pars, sep = "", collapse = ", "), ")", sep = "")))
  }
  
  rescrit <- replicate(M, {
    data <- law.H0(nmax)
    mycrit <- vector("list", nbn)
    for (l in 1:nbstats) {
      if ((length(parstats[[l]]) > 0) && !(any(is.na(parstats[[l]])))) {
        eval(parse(text = paste("mystat <- function(x) vecstats[[l]](x, ", paste(parstats[[l]], collapse = ','), ")", sep = "")))
      } else {
        mystat <- function(x) vecstats[[l]](x)
      }
      for (n in 1:nbn) {
        mycrit[[n]] <- c(mycrit[[n]], mystat(data[1:vectn[n]]))
      }
    }
    mycrit 
  }, simplify = TRUE)

  dim(rescrit) <- c(nbn, M)
  # rescrit is a "matrix" of size nbn x M
  # rescrit[, k] is a list of length nbn
  # rescrit[j, k][[1]] contains nbstats values
  
  quantfunc <- function(x, alter, level) {
    if (alter == 0) ret <- as.vector(quantile(x, probs = c(level / 2, 1 - level / 2)))
    if ((alter == 1) || (alter == 4)) ret <- as.vector(quantile(x, probs = level))
    if ((alter == 2) || (alter == 3)) ret <- as.vector(quantile(x, probs = 1 - level))
    return(ret)
  }
  
  res <- vector("list", nbn)
  for (n in 1:nbn) {
    res[[n]] <- vector("list", nbstats)
    # row s of tmp contains M values of statistic s
    tmp <- matrix(unlist(rescrit[n, ]), nrow = nbstats, ncol = M, byrow = FALSE)
    for (s in 1:nbstats) {
      res[[n]][[s]] <- vector("list", nblevels)
      for (l in 1:nblevels) {
        res[[n]][[s]][[l]] <- quantfunc(tmp[s, ], alter[[s]], levels[l])
      }
    }
  }


  toprint <- matrix(NA, nrow = nbn * nblevels, ncol = 2 + sum(as.vector((unlist(alter) == 0) + 1)))
  for (l in 1:nblevels) {
    for (n in 1:nbn) {
      tmp <- c()
      for (s in 1:nbstats) {
        tmp <- c(tmp, res[[n]][[s]][[l]])
      }
      toprint[n + (l - 1) * nbn, ]  <- c(vectn[n], levels[l], tmp)
    }
  }
  
  return(list(crit = res, toprint = toprint))
}

# Function in pure R to compute the power
powcomp.fastR <- function(vectlaws, stat.indices, vectn = c(20, 50, 100), M = 10 ^ 3, 
                      levels, critval, alter, parstats) {

  eval(parse(text = paste("vecstats <- c(", paste("stat", stat.indices, sep = "", collapse = ", "), ")", sep = "")))

  for (i in 1:length(parstats)) {
    if (any(is.na(parstats[[i]]))) parstats[[i]] <- stat.cstr(as.numeric(strsplit(names(parstats[i]), "stat")[[1]][2]))$stat.pars
  }

  nmax <- max(vectn)
  nbstats <- length(vecstats)
  nbn <- length(vectn)
  nblevels <- length(levels)
  nblaws <- length(vectlaws)
  
  decision <- function(x, alter, crit) {
    if (alter == 0) ret <- mean((x < crit[1]) | (x > crit[2]))
    if ((alter == 1) || (alter == 4)) ret <- mean(x < crit)
    if ((alter == 2) || (alter == 3)) ret <- mean(x > crit)
    return(ret)
  }
  
  res <- vector("list", length(vectlaws))
  for (lawH1 in 1:length(vectlaws)) { 
    res[[lawH1]] <- vector("list", nbn)

    respuis <- replicate(M, {
      data <- vectlaws[[lawH1]](nmax)
      mypuis <- vector("list", nbn)
      for (l in 1:nbstats) {
        if ((length(parstats[[l]]) > 0) && !(any(is.na(parstats[[l]])))) {
          eval(parse(text = paste("mystat <- function(x) vecstats[[l]](x, ", paste(parstats[[l]], collapse = ','), ")", sep = "")))
        } else {
          mystat <- function(x) vecstats[[l]](x)
        }
        for (n in 1:nbn) {
          mypuis[[n]] <- c(mypuis[[n]], mystat(data[1:vectn[n]]))
        }
      }
      mypuis
    }, simplify = TRUE)

    dim(respuis) <- c(nbn, M)
       
    for (n in 1:nbn) {
      res[[lawH1]][[n]] <- vector("list", nbstats)
      tmp <- matrix(unlist(respuis[n, ]), nrow = nbstats, ncol = M, byrow = FALSE)
      for (s in 1:nbstats) {
        res[[lawH1]][[n]][[s]] <- vector("list", nblevels)
        for (l in 1:nblevels) {
          res[[lawH1]][[n]][[s]][[l]] <- decision(tmp[s, ], alter[[s]], critval[[n]][[s]][[l]])
        }
      }
    }
    
    
  }

  toprint <- matrix(NA, nrow = nblaws * nbn * nblevels, ncol = 3 + nbstats)
  for (l in 1:nblevels) {
    for (lawH1 in 1:nblaws) {
      for (n in 1:nbn) {
        tmp <- c()
        for (s in 1:nbstats) {
          tmp <- c(tmp, res[[lawH1]][[n]][[s]][[l]])
        }
        toprint[n + (l - 1) * nblaws * nbn + (lawH1 - 1) * nbn, ]  <- c(lawH1, vectn[n], levels[l], tmp * 100)
      }
    }
  }
  statnames <- c()
  for (stat in stat.indices) statnames <- c(statnames, stat.cstr(stat)$name)
  colnames(toprint) <- c("law", "n", "level", statnames)
  
  return(list(power = res, toprint = toprint))
}


