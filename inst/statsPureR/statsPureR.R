stat2 <- function(x){ # Anderson-Darling; see ad.test in package nortest
    x <- sort(x)
    n <- length(x)
    logp1 <- pnorm((x - mean(x)) / sd(x), log.p = TRUE)
    logp2 <- pnorm(-(x - mean(x)) / sd(x), log.p = TRUE)
    h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
    A <- -n - mean(h)
    AA <- (1 + 0.75 / n + 2.25 / n ^ 2) * A
    if (AA < 0.2) {
        pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA ^ 2)
    }
    else if (AA < 0.34) {
        pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA ^ 2)
    }
    else if (AA < 0.6) {
        pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA ^ 2)
    }
    else if (AA < 10) {
        pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA ^ 2)
    }
    else pval <- 3.7e-24
    res <- list(stat = A, p.value = pval) 
    return(res)
}

stat3 <- function(x) { # Zhang Zc
  n <- length(x)
  phi <- pnorm(sort((x - mean(x)) / sd(x)))
  res <- sum((log((1 / phi - 1) / ((n - 0.5) / ((1:n) - 0.75) - 1))) ^ 2)
  return(res)
}

stat4 <- function(x) { # Zhang Za
  n <- length(x)
  x <- sort(x)
  logp1 <- pnorm((x - mean(x)) / sd(x), lower.tail = TRUE, log.p = TRUE)
  logp2 <- pnorm((x - mean(x)) / sd(x), lower.tail = FALSE, log.p = TRUE)
  res <- -sum(logp1 / (n - (1:n) + 0.5) + logp2 / ((1:n) - 0.5))
  return(10.0 * res - 32.0) # See page 711 in their paper
}

stat5 <- function(x) { # Glen-Leemis-Barr P_s
  n <- length(x)
  Phiz <- pnorm((x - mean(x)) / sd(x))
  Phiz <- sort(Phiz)
  logp1 <- logp2 <- rep(NA, n)
  for (i in 1:n) {
    logp1[i] <- pbeta(Phiz[i], i, n - i + 1, lower.tail = TRUE, log.p = TRUE)
    logp2[i] <- pbeta(Phiz[i], i, n - i + 1, lower.tail = FALSE, log.p = TRUE)
  }
  res <- -n - mean((2 * n + 1 - 2 * (1:n)) * sort(logp1) + (2 * (n:1) - 1) * sort(logp2))
  return(res)
}

stat6 <- function(x) { # D'Agostino-Pearson. Credits to: fBasics package, dagoTest()
    n <- length(x)
    meanX <- mean(x)
    s <- sqrt(mean((x - meanX) ^ 2))
    a3 <- mean((x - meanX) ^ 3)/s ^ 3
    a4 <- mean((x - meanX) ^ 4) / s ^ 4
    SD3 <- sqrt(6 * (n - 2) / ((n + 1) * (n + 3)))
    SD4 <- sqrt(24 * (n - 2) * (n - 3) * n / ((n + 1) ^ 2 * (n + 3) * 
                                              (n + 5)))
    U3 <- a3 / SD3
    U4 <- (a4 - 3 + 6 / (n + 1)) / SD4
    b <- (3 * (n ^ 2 + 27 * n - 70) * (n + 1) * (n + 3)) / ((n - 2) * 
                                                            (n + 5) * (n + 7) * (n + 9))
    W2 <- sqrt(2 * (b - 1)) - 1
    delta <- 1 / sqrt(log(sqrt(W2)))
    a <- sqrt(2 / (W2 - 1))
    Z3 <- delta * log((U3 / a) + sqrt((U3 / a) ^ 2 + 1))
    B <- (6 * (n * n - 5 * n + 2) / ((n + 7) * (n + 9))) * sqrt((6 * 
                                                                 (n + 3) * (n + 5)) / (n * (n - 2) * (n - 3)))
    A <- 6 + (8 / B) * ((2 / B) + sqrt(1 + 4 / (B ^ 2)))
    jm <- sqrt(2 / (9 * A))
    pos <- ((1 - 2 / A) / (1 + U4 * sqrt(2 / (A - 4)))) ^ (1 / 3)
    Z4 <- (1 - 2 / (9 * A) - pos) / jm
    omni <- Z3 ^ 2 + Z4 ^ 2
    pomni <- 1 - pchisq(omni, 2)
    res <- list(stat = omni, p.value = pomni) 
    return(res)
}


stat7 <- function(x){ # Jarque-Bera. Credits to: tseries package, jarque.bera.test()
    n <- length(x)
    m1 <- sum(x) / n
    m2 <- sum((x - m1) ^ 2) / n
    m3 <- sum((x - m1) ^ 3) / n
    m4 <- sum((x - m1) ^ 4) / n
    b1 <- (m3 / m2 ^ (3 / 2)) ^ 2
    b2 <- (m4 / m2 ^ 2)
    STATISTIC <- n * b1 / 6 + n * (b2 - 3) ^ 2 / 24
    PVAL <- 1 - pchisq(STATISTIC, df = 2)
    res <- list(stat = STATISTIC, p.value = PVAL) 
    return(res)
}


stat8 <- function(x) { # Doornik-Hansen. Credits to: normwhn.test package, normality.test1()
    n <- length(x)
    m1 <- sum(x) / n
    v.matrix <- 1 / sqrt(cov.wt(as.matrix(x))$cov)
    xhatprime.matrix <- t(x - m1)
    trprime <- t(v.matrix %*% xhatprime.matrix)
    v1 <- mean(trprime)
    v2 <- (n ^ (-1)) * sum((trprime - v1) ^ 2)
    v3 <- (n ^ (-1)) * sum((trprime - v1) ^ 3)
    v4 <- (n ^ (-1)) * sum((trprime - v1) ^ 4)
    rtb1 <- v3 / (v2 ^ (3 / 2))
    b2 <- v4 / v2 ^ 2
    beta <- (3 * (n ^ 2 + 27 * n - 70) * (n + 1) * (n + 3)) / ((n - 
        2) * (n + 5) * (n + 7) * (n + 9))
    w2 <- (-1) + sqrt(2 * (beta - 1))
    delta <- 1 / sqrt(log(sqrt(w2)))
    f <- (w2 - 1) / 2
    g <- (n + 1) * (n + 3) / (6 * (n - 2))
    h <- sqrt(f * g)
    y <- rtb1 * h
    z1 <- delta * log(y + sqrt(y ^ 2 + 1))
    del <- ((n - 3) * (n + 1) * (n ^ 2 + 15 * n - 4))
    aye <- ((n - 2) * (n + 5) * (n + 7) * (n ^ 2 + 27 * n - 70)) / (6 * 
        del)
    cee <- ((n - 7) * (n + 5) * (n + 7) * (n ^ 2 + 2 * n - 5)) / (6 * 
        del)
    alp <- aye + ((rtb1 ^ 2) * cee)
    kap <- ((n + 5) * (n + 7) * (n ^ 3 + 37 * n ^ 2 + 11 * n - 313)) / (12 * 
        del)
    chi <- (b2 - 1 - rtb1 ^ 2) * (2 * kap)
    chi <- abs(chi)
    z2 <- (((chi / (2 * alp)) ^ (1 / 3)) - 1 + (1 / ((9 * alp)))) * sqrt(9 * 
        alp)
    stat <- z1 ^ 2 + z2 ^ 2
    pval <- 1 - pchisq(stat, 2)
    res <- list(stat = stat, p.value = pval) 
    return(res)
}


stat17 <-function (x) { # Bonett-Seier. Credits to moments package, bonett.test()
    x <- sort(x)
    n <- length(x)
    rho <- sqrt(sum((x - mean(x)) ^ 2) / n)
    tau <- sum(abs(x - mean(x))) / n
    omega <- 13.29 * (log(rho) - log(tau))
    z <- sqrt(n + 2) * (omega - 3) / 3.54
    pval <- 2* pnorm(z, lower.tail = FALSE)
        if (pval > 1) 
            pval <- 2 - pval
    res <- list(stat = z, p.value = pval) 
    return(res)
}


stat26 <- function(x){ # Chen-Shapiro
  n <- length(x)
  m <- qnorm(((1:n) - 0.375) / (n + 0.25))
  xo <- sort(x)
  QH <- sum(diff(xo) / diff(m)) / (n - 1) / sd(x); 
  res <- sqrt(n) * (1 - QH)
  return(res)
}

stat29 <- function(x){ # BCMR
  n <- length(x)
  xo <- sort(x)
  dt <- 0.000001
  k <- 1:n
  sigma <- 0
  for (i in k){
    u <- seq((k[i] - 1) / n + dt / 2, k[i] / n - dt / 2, dt)
    sigma <- sigma +  sum(qnorm(u)) * xo[i] * dt
  }
  stat <- 1 - sigma ^ 2 / (var(x) * (n - 1) / n)
  return(stat)
}

stat30 <- function(x) { # Coin. Credits to http://digilander.libero.it/polynomial/B3test.R
    normorder <- function(n){
        dt <- 0.01
        t <- seq(-10, 10, dt)
        int <- function(r) {
            rinv <- n - r + 1
            logh <- lgamma(n + 1) - lgamma(rinv) - lgamma(n - rinv + 1) +
                (rinv - 1) * pnorm(t, low = FALSE, log.p = TRUE) + (n - rinv) * pnorm(t, low = TRUE, log.p = TRUE) +
                dnorm(t, log = TRUE) 
            h <- t * exp(logh)
            return(sum(h) * dt)
        }
        res <- apply(cbind(1:n), 1, int)
        return(res)
    }
    n <- length(x)
    a <- normorder(n)
    x <- sort(x)
    z <- (x - mean(x)) / sd(x)
    mod <- lm(z ~ 0 + a + I(a ^ 3))
    stat <- mod$coef[2] ^ 2
    return(stat = stat)
}

stat35 <- function(x) { # GLB^*
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  uhat <- 0.5 * exp(-abs(x - muhat) / sighat)
  uhat[x >= muhat] <- 1 - uhat[x >= muhat]
  pihat <- rep(NA, n)
  for (i in 1:n) {
    pihat[i] <- pbeta(uhat[i], i, n - i + 1)
  }
  pihat <- sort(pihat)
  res <- - n - mean((2 * n + 1 - 2 * (1:n)) * log(pihat) + (2 * (1:n) - 1) * log(1 - pihat))
  return(res)
}

stat36 <- function(x) { # BRT_3
  n <- length(x)
  xbar <- mean(x)
  sn2 <- mean((x - xbar) ^ 2)
  h3 <- function(z) (-12 * z + z ^ 3) / (12 * sqrt(3))
  res <- sum(h3((x - xbar) / sqrt(0.5 * sn2))) / sqrt(n)
  return(res)
}

stat37 <- function(x) { # BRT_4
  n <- length(x)
  xbar <- mean(x)
  sn2 <- mean((x - xbar) ^ 2)
  h4 <- function(z) (43.2 -33.6 * z ^ 2 + z ^ 4) / (24 * sqrt(149 / 5))
  res <- sum(h4((x - xbar) / sqrt(0.5 * sn2))) / sqrt(n)
  return(res)
}

stat38 <- function(x) { # Zepd
    n <- length(x) 
    euler <- -digamma(1)
    z <- (x - mean(x)) / sqrt(var(x) * (n - 1) / n) 
    z <- z[! (z == 0)] 
    K2 <- sum(z ^ 2 * log(abs(z))) / n 
    alpha.n <- -0.06 + 2.1 / n ^ 0.67
    K2.box.cox <- ((2.0 * K2) ^ alpha.n - 1.0) / alpha.n
    esp.K2.box.cox <- - ((2.0 - log(2.0) - euler) ^ (-0.06) - 1.0) / 0.06 -
        1.32 / n ^ 0.95
    var.K2.box.cox <- (1.0 / n) * ((2.0 - log(2.0) - euler) ^ (-2.12) *
                                   (3.0 * pi ^ 2 - 28.0) / 2.0 - 3.78 / n ^ 0.733)
    stat.Z.K2 <- (K2.box.cox - esp.K2.box.cox) / sqrt(var.K2.box.cox)
    p.value <- 2.0 * pnorm(abs(stat.Z.K2), lower.tail = FALSE)
    res <- list(test = "2nd-power kurtosis directional normality test Zepd",
                stat.Z.K2 = stat.Z.K2, p.value = p.value)
    return(res)
}

stat39 <- function(x) { # Xapd
    n <- length(x) 
    euler <- - digamma(1)
    z <- (x - mean(x)) / sqrt(var(x) * (n - 1) / n) 
    z <- z[! (z == 0)] 
    B2 <- sum(z ^ 2 * sign(z)) / n 
    K2 <- sum(z ^ 2 * log(abs(z))) / n 
    var.B2 <- (1 / n) * (3.0 - 8.0 / pi) * (1.0 - 1.9 / n) 
    Z.B2 <-  B2 / sqrt(var.B2)
    esp.net.K2 <- ((2.0 - log(2.0) - euler) / 2.0) ^ (1.0 / 3.0) *
        (1 - 1.026 / n) 
    var.net.K2 <- (1.0 / n) * (1.0 / 72.0) *
        ((2.0 - log(2.0) - euler) / 2.0) ^ (-4.0 / 3.0) *
        (3.0 * pi ^ 2 - 28.0) * (1.0 - 2.25 / n ^ 0.8) ;
    Z.net.K2 <- ((K2 - B2 ^ 2) ^ (1.0 / 3.0) - esp.net.K2) /
        sqrt(var.net.K2)
    Xapd.stat <- Z.B2 ^ 2 + Z.net.K2 ^ 2 ;
    p.value <- pchisq(Xapd.stat, 2, lower.tail = FALSE) ;
    res <- list(test = "2nd-power skewness and kurtosis omnibus normality test Xapd",
                B2 = B2, K2 = K2, Z.B2 = Z.B2, Z.net.K2 = Z.net.K2,
                Xapd.stat = Xapd.stat, p.value = p.value)
    return(res) 
}

stat42 <- function(x, theta1, theta2, lambda) { # AD
  n <- length(x)
  islambdaone <- abs(lambda - 1.0) < 10 ^ -15
  islambdatwo <- abs(lambda - 2.0) < 10 ^ -15
  if (islambdaone) { # lambda = 1
    muhat <- median(x)
  } else if (islambdatwo) { # lambda = 2
    muhat <- mean(x)
  } else {
    f.opt <- function(mu, lambda, x) {
      sum(abs(x - mu) ^ (lambda - 1.0) * sign(x - mu))
    }
    muhat <- uniroot(f = f.opt, interval = c(min(x), max(x)),
                     lambda = lambda, x = x, 
                     tol = 10 ^ -12, #.Machine$double.eps ^ 0.25,
                     maxiter = 1000, trace = 0)$root
  }
  sighat <- (mean(abs(x - muhat) ^ lambda)) ^ (1.0 / lambda)
  pAPD <- function(q, theta1, theta2, lambda, mu = 0, sigma = 1){
    z <- (q - mu) / sigma
    delta <- 2 * theta1 ^ theta2 * (1 - theta1) ^ theta2 / (theta1 ^ theta2 + (1 - theta1) ^ theta2)
    F <- theta1 * (1 - pgamma(delta / lambda * (-pmin(0, z) / theta1) ^ theta2, shape = 1 / theta2)) +
         (1 - theta1) * pgamma(delta / lambda * (pmax(0, z) / (1 - theta1)) ^ theta2, shape = 1 / theta2)
    return(F) 
  }
  uhat <- pAPD(x, theta1, theta2, lambda, muhat, sighat)
  res <- -n - mean((2 * (1:n) - 1) * log(uhat) + (2 * (n:1) - 1) * log(1 - uhat))
  return(res)
}

stat43 <- function(x, theta1, theta2, lambda) { # CvM
  n <- length(x)
  islambdaone <- abs(lambda - 1.0) < 10 ^ -15
  islambdatwo <- abs(lambda - 2.0) < 10 ^ -15
  if (islambdaone) { # lambda = 1
    muhat <- median(x)
  } else if (islambdatwo) { # lambda = 2
    muhat <- mean(x)
  } else {
    f.opt <- function(mu, lambda, x) {
      sum(abs(x - mu) ^ (lambda - 1.0) * sign(x - mu))
    }
    muhat <- uniroot(f = f.opt, interval = c(min(x), max(x)),
                     lambda = lambda, x = x, 
                     tol = 10 ^ -12, #.Machine$double.eps ^ 0.25,
                     maxiter = 1000, trace = 0)$root
  }
  sighat <- (mean(abs(x - muhat) ^ lambda)) ^ (1.0 / lambda)
  pAPD <- function(q, theta1, theta2, lambda, mu = 0, sigma = 1){
    z <- (q - mu) / sigma
    delta <- 2 * theta1 ^ theta2 * (1 - theta1) ^ theta2 / (theta1 ^ theta2 + (1 - theta1) ^ theta2)
    F <- theta1 * (1 - pgamma(delta / lambda * (-pmin(0, z) / theta1) ^ theta2, shape = 1 / theta2)) +
         (1 - theta1) * pgamma(delta / lambda * (pmax(0, z) / (1 - theta1)) ^ theta2, shape = 1 / theta2)
    return(F) 
  }
  uhat <- pAPD(x, theta1, theta2, lambda, muhat, sighat)
  res <- 1 / (12 * n) + sum(((2 * (1:n) - 1) / (2 * n) - uhat) ^ 2)
  return(res)
}

stat44 <- function(x, theta1, theta2, lambda) { # Wa
  n <- length(x)
  islambdaone <- abs(lambda - 1.0) < 10 ^ -15
  islambdatwo <- abs(lambda - 2.0) < 10 ^ -15
  if (islambdaone) { # lambda = 1
    muhat <- median(x)
  } else if (islambdatwo) { # lambda = 2
    muhat <- mean(x)
  } else {
    f.opt <- function(mu, lambda, x) {
      sum(abs(x - mu) ^ (lambda - 1.0) * sign(x - mu))
    }
    muhat <- uniroot(f = f.opt, interval = c(min(x), max(x)),
                     lambda = lambda, x = x, 
                     tol = 10 ^ -12, #.Machine$double.eps ^ 0.25,
                     maxiter = 1000, trace = 0)$root
  }
  sighat <- (mean(abs(x - muhat) ^ lambda)) ^ (1.0 / lambda)
  pAPD <- function(q, theta1, theta2, lambda, mu = 0, sigma = 1){
    z <- (q - mu) / sigma
    delta <- 2 * theta1 ^ theta2 * (1 - theta1) ^ theta2 / (theta1 ^ theta2 + (1 - theta1) ^ theta2)
    F <- theta1 * (1 - pgamma(delta / lambda * (-pmin(0, z) / theta1) ^ theta2, shape = 1 / theta2)) +
         (1 - theta1) * pgamma(delta / lambda * (pmax(0, z) / (1 - theta1)) ^ theta2, shape = 1 / theta2)
    return(F) 
  }
  uhat <- pAPD(x, theta1, theta2, lambda, muhat, sighat)
  CvMstar <- 1 / (12 * n) + sum(((2 * (1:n) - 1) / (2 * n) - uhat) ^ 2)
  res <- CvMstar - n * (mean(uhat) - 0.5) ^ 2
  return(res)
}

stat45 <- function(x, theta1, theta2, lambda) { # \sqrt{n}D = KS
  n <- length(x)
  islambdaone <- abs(lambda - 1.0) < 10 ^ -15
  islambdatwo <- abs(lambda - 2.0) < 10 ^ -15
  if (islambdaone) { # lambda = 1
    muhat <- median(x)
  } else if (islambdatwo) { # lambda = 2
    muhat <- mean(x)
  } else {
    f.opt <- function(mu, lambda, x) {
      sum(abs(x - mu) ^ (lambda - 1.0) * sign(x - mu))
    }
    muhat <- uniroot(f = f.opt, interval = c(min(x), max(x)),
                     lambda = lambda, x = x, 
                     tol = 10 ^ -12, #.Machine$double.eps ^ 0.25,
                     maxiter = 1000, trace = 0)$root
  }
  sighat <- (mean(abs(x - muhat) ^ lambda)) ^ (1.0 / lambda)
  pAPD <- function(q, theta1, theta2, lambda, mu = 0, sigma = 1){
    z <- (q - mu) / sigma
    delta <- 2 * theta1 ^ theta2 * (1 - theta1) ^ theta2 / (theta1 ^ theta2 + (1 - theta1) ^ theta2)
    F <- theta1 * (1 - pgamma(delta / lambda * (-pmin(0, z) / theta1) ^ theta2, shape = 1 / theta2)) +
         (1 - theta1) * pgamma(delta / lambda * (pmax(0, z) / (1 - theta1)) ^ theta2, shape = 1 / theta2)
    return(F) 
  }
  uhat <- pAPD(x, theta1, theta2, lambda, muhat, sighat)
  res <- sqrt(n) * max(max(uhat - ((1:n) - 1) / n), max((1:n) / n - uhat))
  return(res)
}

stat46 <- function(x, theta1, theta2, lambda) { # \sqrt{n}V = Ku
  n <- length(x)
  islambdaone <- abs(lambda - 1.0) < 10 ^ -15
  islambdatwo <- abs(lambda - 2.0) < 10 ^ -15
  if (islambdaone) { # lambda = 1
    muhat <- median(x)
  } else if (islambdatwo) { # lambda = 2
    muhat <- mean(x)
  } else {
    f.opt <- function(mu, lambda, x) {
      sum(abs(x - mu) ^ (lambda - 1.0) * sign(x - mu))
    }
    muhat <- uniroot(f = f.opt, interval = c(min(x), max(x)),
                     lambda = lambda, x = x, 
                     tol = 10 ^ -12, #.Machine$double.eps ^ 0.25,
                     maxiter = 1000, trace = 0)$root
  }
  sighat <- (mean(abs(x - muhat) ^ lambda)) ^ (1.0 / lambda)
  pAPD <- function(q, theta1, theta2, lambda, mu = 0, sigma = 1){
    z <- (q - mu) / sigma
    delta <- 2 * theta1 ^ theta2 * (1 - theta1) ^ theta2 / (theta1 ^ theta2 + (1 - theta1) ^ theta2)
    F <- theta1 * (1 - pgamma(delta / lambda * (-pmin(0, z) / theta1) ^ theta2, shape = 1 / theta2)) +
         (1 - theta1) * pgamma(delta / lambda * (pmax(0, z) / (1 - theta1)) ^ theta2, shape = 1 / theta2)
    return(F) 
  }
  uhat <- pAPD(x, theta1, theta2, lambda, muhat, sighat)
  res <- sqrt(n) * (max(uhat - ((1:n) - 1) / n) + max((1:n) / n - uhat))
  return(res)
}

stat48 <- function(x, a = 2) { # Me_{a}^{(1)}
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  zhat <- (x - muhat) / sighat
  z2 <- zhat ^ 2
  a2 <- a ^ 2
  diff.z2 <- outer(zhat, zhat, FUN = "-") ^ 2
  res <- (2 * n) / a -
    4 * a * sum(1 / (a2 + z2) + 2 * (a2 - 3 * z2) / (a2 + z2) ^ 3) +
    2 * a * sum(1 / (a2 + diff.z2) + 24 * (a2 ^ 2 + 5 * diff.z2 ^ 2 - 10 * a2 * diff.z2) / (a2 + diff.z2) ^ 5 + 4 * (a ^ 2 - 3 * diff.z2) / (a2 + diff.z2) ^ 3) / n
  return(res)
}

stat50 <- function(x, a) { # Me_{a}^{(2)}
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  zhat <- (x - muhat) / sighat
  z2 <- zhat ^ 2
  a2 <- a ^ 2
  diff.z2 <- outer(zhat, zhat, FUN = "-") ^ 2
  pia <- sqrt(pi / a)
  res <- n * pia -
    2 * pia * sum((1 - (z2 - 2 * a) / (4 * a2)) * exp(-z2 / (4 * a))) +
    2 * pia * sum((0.5 + (diff.z2 ^ 2 + 12 * a2 - 12 * a * diff.z2) / (32 * a2 ^ 2) - (diff.z2 - 2 * a) / (4 * a2)) * exp(-diff.z2 / (4 * a))) / n
  return(res)
}

stat51 <- function(x, version = 1, m = NULL) { # CK_v, CK_e, CK_c
  n <- length(x)
  if (is.null(m) || (m  == 0)) {
    if (version == 1) {
      if (n >= 51) {
        m  <-  ceiling(n / 10)
      } else {
        mTV <-  c(1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
                  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0,
                  4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                  5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 6.0, 6.0, 6.0, 6.0)
        m  <-  mTV[n - 3]
      }
    } else if (version == 2) {
      if (n >= 51) {
        m  <-  ceiling(n / 10)
      } else {
        mTC <-  c(1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 3.0, 2.0, 3.0, 3.0, 3.0,
                  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0,
                  4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0,
                  5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0)
        m  <-  mTC[n - 3]
      }
    } else if (version == 3) {
      if (n >= 12) {
        m  <-  2
      } else {
        mTE <-  c(1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0)
        m  <-  mTE[n - 3]
      }
    }
  }
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  min <- x[1]
  max <- x[n]
  CKv <- 0
  for (i in 1:n) {
    if ((i + m) <= n) {temp1 <- x[i + m]} else {temp1 <- max}
    if ((i - m) >= 1) {temp2 <- x[i - m]} else {temp2 <- min}
    CKv  <- CKv + log(temp1 - temp2)
  }
  CKv <- n * exp(CKv / n) / (2 * m * sighat)
  CKe <- 0
  for (i in 1:(n - m)) {
    if ((i + m) <= n) {temp1 <- x[i + m]} else {temp1 <- max}
    CKe  <- CKe + log(temp1 - x[i])
  }
  CKe <- exp(CKe / (n - m)) * exp(sum(1/(m:n))) / sighat
  prod <- 0
  for (i in 1:n) {
    tmp <- 0
    for (k in (i - m):(i + m)) {
      if (k > n) {
        temp <- max
      } else if (k < 1) {
        temp <- min
      } else {
        temp <- x[k]
      }
      tmp <- tmp + temp
    }
    tmp <- tmp / (2 * m + 1)
    temp1 <- temp2 <- 0
    for (j in (i - m):(i + m)) {
      if (j > n) {
        temp <- max
      } else if (j < 1) {
        temp <- min
      } else {
        temp <- x[j]
      }
      temp1 <- temp1 + (j - i) * (temp - tmp) / n
      temp2 <- temp2 + (temp - tmp) ^ 2
    }
    prod  <- prod + log(temp1 / temp2)
  }
  CKc <- exp(-prod / n) / sighat
  if (version == 1) res <- CKv
  if (version == 2) res <- CKc
  if (version == 3) res <- CKe
  return(res)
}

stat55 <- function(x) { # V_3 = RB_3
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  h3 <- function(z) (-12 * z + z ^ 3) / (12 * sqrt(3))
  res <- sum(h3((x - muhat) / sighat)) / sqrt(n)
  return(res)
}


stat56 <- function(x) { # V_4 = RB_4
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  h4 <- function(z) (43.2 -33.6 * z ^ 2 + z ^ 4) / (24 * sqrt(149 / 5))
  res <- sum(h4((x - muhat) / sighat)) / sqrt(n)
  return(res)
}



stat57 <- function(x) { # LK
  n <- length(x)
  xbar <- mean(x)
  sn <- sqrt(mean((x - xbar) ^ 2))
  Phi <- function(z) {
    rez <- 0.5 * exp(-abs(z))
    rez[z >= 0] <- 1 - rez[z >= 0]
    return(rez)
  }
  res <- 1.856 * n * ((mean(cos(2 * pi * Phi(sqrt(2) * (x - xbar) / sn)))) ^ 2 + (mean(sin(2 * pi * Phi(sqrt(2) * (x - xbar) / sn)))) ^ 2)
  return(res)
}

stat59 <- function(x) { # BS
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  y <- sort(abs(x - muhat) / sighat)
  w <- (n - (1:n) + 1) * (y - c(0, y[-n]))
  v <- cumsum(w)[-n] / sum(w)
  vbar <- mean(v)
  res <- 12 * (n - 1) * (vbar - 0.5) ^ 2 + 5 * (n - 1) * (n - 2 + 6 * n * vbar - 12 * sum((1:(n - 1)) * v) / (n - 1)) ^ 2 / ((n + 2) * (n - 2))
  return(res)
}

stat60 <- function(x) { # Ge
  n <- length(x)
  xbar <- mean(x)
  sighat <- mean(abs(x - median(x)))
  sqrtv1 <- mean((x - xbar) ^ 3) / (2 * sqrt(2) * sighat ^ 3)
  v2 <- mean((x - xbar) ^ 4) / (4 * sighat ^ 4)
  C1 <- 60
  C2 <- 1200
  res <- n * (sqrtv1 ^ 2 / C1 + (v2 - 6) ^ 2 / C2)
  return(res)
}

stat84 <- function(x, group, version, m = NULL) { # AP_{l} and AP_{l}^{(MLE)} for version = l = v, e, y, a, z
  n <- length(x)
  if (is.null(m) || (m == 0)) {
    if (group == 1 && version == 5) {
      if ((1 <= n) && (n <= 8)) m <- 1
      if ((9 <= n) && (n <= 25)) m <- 2
      if ((26 <= n) && (n <= 40)) m <- 3
      if ((41 <= n) && (n <= 60)) m <- 4
      if ((61 <= n) && (n <= 90)) m <- 5
      if ((91 <= n) && (n <= 120)) m <- 6
      if (121 <= n) m <- round((log(n)) ^ 1.15)
    } else {
      if ((1 <= n) && (n <= 8)) m <- 1
      if ((9 <= n) && (n <= 15)) m <- 2
      if ((16 <= n) && (n <= 25)) m <- 4
      if ((26 <= n) && (n <= 40)) m <- 5
      if ((41 <= n) && (n <= 60)) m <- 6
      if ((61 <= n) && (n <= 90)) m <- 7
      if ((91 <= n) && (n <= 120)) m <- 8
      if (121 <= n) m <- round((log(n)) ^ 1.35)
    }
  }
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  min <- x[1]
  max <- x[n]
  Hv <- 0
  for (i in 1:n) {
    if ((i + m) <= n) {temp1 <- x[i + m]} else {temp1 <- max}
    if ((i - m) >= 1) {temp2 <- x[i - m]} else {temp2 <- min}
    Hv = Hv + log(n * (temp1 - temp2) / (2 * m))
  }
  Hv <- Hv / n
  xi <- rep(NA, n + 1)
  for (i in 1:m) {
    xi[i] <- (m + 2 - i) * min / (2 * m) + sum(x[2:(i + m - 1)]) / (2 * m)
  }
  for (i in (m + 1):(n - m + 1)) {
    xi[i] <- sum(x[(i - m):(i + m - 1)]) / (2 * m)
  }
  for (i in (n - m + 2):(n + 1)) {
    xi[i] <- (i - n + m) * max / (2 * m) + sum(x[(i - m):(n - 1)]) / (2 * m)
  }
  ci <- rep(NA, n)
  for (i in 1:m) {
    ci[i] <- 1 + (i - 1) / m
  }
  for (i in (m+ 1):(n - m)) {
    ci[i] <- 2
  }
  for (i in (n - m + 1):n) {
    ci[i] <- 1 + (n - i) / m
  }
  He <- 0
  for (j in 1:n) {
    if ((j + m) <= n) {temp1 <- x[j + m]} else {temp1 <- max}
    if ((j - m) >= 1) {temp2 <- x[j - m]} else {temp2 <- min}
    He = He + log(n * (temp1 - temp2) / (ci[j] * m))
  }
  He <- He / n
  eta <- rep(NA, n + 1)
  for (i in 1:m) {
    eta[i] <- xi[m + 1] - sum((x[(m + i):(2 * m)] - min) / (m + (i:m) - 1))
  }
  for (i in (m + 1):(n - m + 1)) {
    eta[i] <- xi[i]
  }
  for (i in (n - m + 2):(n + 1)) {
    eta[i] <- xi[n - m + 1] + sum((max - x[(n - 2 * m + 1):(i - m - 1)]) / (n + m - ((n - m + 2):i) + 1))
  }
  Fhat <- function(i) {
    if (i == 1) rez <- n / (n - 1) + if (abs(x[1] - x[2]) < 10 ^ -50) 0 else n / (2 * n - 1)
    if ((i >= 2) && (i <= (n - 1))) rez <- (i * (n - 1) + 1) / (n - 1) + if (abs(x[i + 1] - x[i - 1]) < 10 ^ -50) 0 else (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
    if (i == n) rez <- (n ^ 2 - 2 * n + 2) / (n - 1) + if (abs(x[n - 2] - x[n]) < 10 ^ -50) 0 else if ((abs(x[n - 1] - x[n]) < 10 ^ -50) && (abs(x[n - 2] - x[n - 1]) > 10 ^ -50)) 1 else 1 + (n - 1) / (2 * n - 1)
    return((n - 1) * rez / (n * (n + 1)))
  }
  sumFhat <- 0
  for (j in 1:n) {
    if ((j + m) <= n) {i1 <- j + m} else {i1 <- n}
    if ((j - m) >= 1) {i2 <- j - m} else {i2 <- 1}
    sumFhat  <-  sumFhat + Fhat(i1) - Fhat(i2)
  }
  Hy <- 0
  for (j in 1:n) {
    if ((j + m) <= n) {temp1 <- x[j + m]; i1 <- j + m} else {temp1 <- max; i1 <- n}
    if ((j - m) >= 1) {temp2 <- x[j - m]; i2 <- j - m} else {temp2 <- min; i2 <- 1}
    Hy = Hy + ((Fhat(i1) - Fhat(i2)) / sumFhat) * log((temp1 - temp2) / (Fhat(i1) - Fhat(i2)))
  }
  zeta <- rep(NA, n + 1)
  for (i in 1:(n + 1)) {
    zeta[i] <- 2 * m * xi[i] / sumFhat
  }
  ai <- rep(1, n)
  for (i in (m + 1):(n - m)) {
    ai[i] <- 2
  }
  Ha <- 0
  for (j in 1:n) {
    if ((j + m) <= n) {temp1 <- x[j + m]} else {temp1 <- max}
    if ((j - m) >= 1) {temp2 <- x[j - m]} else {temp2 <- min}
    Ha = Ha + log(n * (temp1 - temp2) / (ai[j] * m))
  }
  Ha <- Ha / n
  nu <- rep(NA, n + 1)
  for (i in 1:m) {
    nu[i] <- xi[m + 1] - sum(x[(m + i):(2 * m)] - min) / m
  }
  for (i in (m + 1):(n - m + 1)) {
    nu[i] <- xi[i]
  }
  for (i in (n - m + 2):(n + 1)) {
    nu[i] <- xi[n - m + 1] + sum(max - x[(n - 2 * m + 1):(i - m - 1)]) / m
  }
  zi <- rep(NA, n)
  for (i in 1:m) {
    zi[i] <- i / m
  }
  for (i in (m + 1):(n - m)) {
    zi[i] <- 2
  }
  for (i in (n - m + 1):n) {
    zi[i] <- (n - i + 1) / m
  }
  Hz <- 0
  for (j in 1:n) {
    if ((j + m) <= n) {temp1 <- x[j + m]} else {temp1 <- max}
    if ((j - m) >= 1) {temp2 <- x[j - m]} else {temp2 <- min}
    Hz = Hz + log(n * (temp1 - temp2) / (zi[j] * m))
  }
  Hz <- Hz / n
  tau <- rep(NA, n + 1)
  for (i in 1:m) {
    tau[i] <- xi[m + 1] - sum((x[(m + i):(2 * m)] - min) / (i:m))
  }
  for (i in (m + 1):(n - m + 1)) {
    tau[i] <- xi[i]
  }
  for (i in (n - m + 2):(n + 1)) {
    tau[i] <- xi[n - m + 1] + sum((max - x[(n - 2 * m + 1):(i - m - 1)]) / (n - (n - m + 2):i + 2))
  }
  thetahat <- function(beta) {
    if (n %% 2 == 0) { # n even
      rez <- 0
      for (i in 1:(n / 2)) {
        rez <- rez + 0.5 * (beta[i] + beta[i + 1])
      }
      rez <- -rez / n
      for (i in ((n / 2) + 1):n) {
        rez <- rez + 0.5 * (beta[i] + beta[i + 1]) / n
      }
    } else { # n odd
      term1 <- 0
      for (i in 1:((n - 1) / 2)) {
        term1 <- term1 + 0.5 * (beta[i] + beta[i + 1])
      }
      term2 <- 0.25 * (beta[(n + 1) / 2 + 1] - beta[(n + 1) / 2])
      term3 <- 0
      for (i in ((n + 1) / 2 + 1):n) {
        term3 <- term3 + 0.5 * (beta[i] + beta[i + 1])
      }
      rez <- (-term1 + term2 + term3) / n
    }
    return(rez)
  }
  APv <- log(2 * thetahat(xi)) + 1 - Hv
  APe <- log(2 * thetahat(eta)) + 1 - He
  APy <- log(2 * thetahat(zeta)) + 1 - Hy
  APa <- log(2 * thetahat(nu)) + 1 - Ha
  APz <- log(2 * thetahat(tau)) + 1 - Hz
  APvMLE <- log(2 * sighat) + 1 - Hv  
  APeMLE <- log(2 * sighat) + 1 - He  
  APyMLE <- log(2 * sighat) + 1 - Hy  
  APaMLE <- log(2 * sighat) + 1 - Ha  
  APzMLE <- log(2 * sighat) + 1 - Hz  

  if (group == 1) {
    if (version == 1) {
      res <- APv
    }
    if (version == 2) {
      res <- APe
    }
    if (version == 3) {
      res <- APy
    }
    if (version == 4) {
      res <- APa
    }
    if (version == 5) {
      res <- APz
    }
  }
  if (group == 2) {
    if (version == 1) {
      res <- APvMLE
    }
    if (version == 2) {
      res <- APeMLE
    }
    if (version == 3) {
      res <- APyMLE
    }
    if (version == 4) {
      res <- APaMLE
    }
    if (version == 5) {
      res <- APzMLE
    }
  }
  return(res)
}

stat91 <- function(x) { # GV_1
  n <- length(x)
  xbar <- mean(x)
  sn2 <- mean((x - xbar) ^ 2)
  res <- sqrt(4 * n) * (sqrt(0.5 * sn2) / mean(abs(x - xbar)) - 1)
  return(res)
}

stat92 <- function(x) { # GV_2
  n <- length(x)
  xbar <- mean(x)
  sn2 <- mean((x - xbar) ^ 2)
  sighat <- mean(abs(x - median(x)))
  res <- sqrt(4 * n) * (sqrt(0.5 * sn2) / sighat - 1)
  return(res)
}

stat93 <- function(x) { # Ho_K
  n <- length(x)
  xbar <- mean(x)
  sn2 <- mean((x - xbar) ^ 2)
  res <- mean((x - xbar) ^ 4) / sn2 ^ 2
  return(res)
}

stat94 <- function(x) { # Ho_U
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  xbar <- mean(x)
  sn2 <- mean((x - xbar) ^ 2)
  res <- sqrt(sn2) / sighat
  return(res)
}

stat95 <- function(x) { # Ho_V
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  res <- 0.5 * (x[n] - x[1]) / sighat
  return(res)
}

stat96 <- function(x) { # Ho_W
  n <- length(x)
  x <- sort(x)
  xbar <- mean(x)
  sn2 <- mean((x - xbar) ^ 2)
  res <- 0.5 * (x[n] - x[1]) / sqrt(sn2)
  return(res)
}

stat97 <- function(x) { # SR^*
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  zhat <- (x - muhat) / sighat
  res <- 2 * sum(abs(zhat) + exp(-abs(zhat))) - 1.5 * n - 2 * mean((2 * (1:n) - 1 - n) * zhat)
  return(res)
}

stat98 <- function(x, version) { # AB_{x} with x = KL, He, Je, TV or \chi^2
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  hstar <- sighat * (2 / (n * sqrt(pi))) ^ (1 / 5)
  K <- function(z) return(exp(-0.5 * z ^ 2) / sqrt(2 * pi))
  fhat <- function(y) {
    rez <- rep(NA, length(y))
    for (i in 1:length(y)) {
      rez[i] <- mean(K((y[i] - y) / hstar)) / hstar
    }
    return(rez)
  }
  f <- function(y) 0.5 * exp(-abs((y - muhat) / sighat)) / sighat
  if (version == 1) {
    nu <- function(t) return(t * log(t))
  } else if (version == 2) {
    nu  <- function(t) 0.5 * (sqrt(t) - 1) ^ 2
  } else if (version == 3) {
    nu  <- function(t) (t - 1) * log(t)
  } else if (version == 4) {
    nu <- function(t) abs(t - 1)
  } else if (version == 5) {
    nu <- function(t) (t - 1) ^ 2
  }
  tmp <- fhat(x) / f(x)
  res <- mean(nu(tmp) / tmp)
  return(res)
}

stat99 <- function(x, delta = 0.5) { # A_{ratio}
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  min <- x[1]
  max <- x[n]
  term2 <- 0.5 * exp(-abs(x - muhat) / sighat) / sighat
  tmp  <- 1 # m = 1 here
  for (j in 1:n) {
    if ((j + 1) <= n) {temp1 <- x[j + 1]} else {temp1 <- max}
    if (j >= 2) {temp2 <- x[j - 1]} else {temp2 <- min}
    tmp  <- tmp * 2 / (n * (temp1 - temp2))
    tmp <- tmp / term2[j]
  }
  term1  <-  tmp
  #cte <- min(as.integer(n ^ delta), as.integer(n / 2)) - 1
  bound <- min(n ^ delta, n / 2)
  cte <- as.integer(bound) - (bound - as.integer(bound) < .Machine$double.eps)
  if (cte >= 2) {
    for (m in 2:cte) {
      tmp  <- 1
      for (j in 1:n) {
        if ((j + m) <= n) {temp1 <- x[j + m]} else {temp1 <- max}
        if ((j - m) >= 1) {temp2 <- x[j - m]} else {temp2 <- min}
        tmp  <- tmp * 2 * m / (n * (temp1 - temp2))
        tmp <- tmp / term2[j]
      }
      if (tmp < term1) term1 <- tmp
    }
  }
  res <- term1
  return(res)
}

stat100 <- function(x, m = NULL) { # A_{entropy}
  n <- length(x)
  if (is.null(m) || (m == 0)) {
    if (n <= 3) m <- 1
    if ((n == 4) || (n == 5)) m <- 2
    if (n >= 6) m <- round((n + 2) / 5)
  }
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  min <- x[1]
  max <- x[n]
  Phi <- function(z) {
    rez <- 0.5 * exp(-abs(z))
    rez[z >= 0] <- 1 - rez[z >= 0]
    return(rez)
  }
  sum <- 0
  for (j in 1:n) {
    if ((j + m) <= n) {temp1 <- x[j + m]} else {temp1 <- max}
    if ((j - m) >= 1) {temp2 <- x[j - m]} else {temp2 <- min}
    sum <- sum + log(n * (Phi((temp1 - muhat) / sighat) - Phi((temp2 - muhat) / sighat)) / (2 * m))
  }
  res <- -sum / n
  return(res)
}

stat101 <- function(x, version, theta1, theta2, lambda) { # Z_K, Z_A, Z_C
  n <- length(x)
  islambdaone <- abs(lambda - 1.0) < 10 ^ -15
  islambdatwo <- abs(lambda - 2.0) < 10 ^ -15
  if (islambdaone) { # lambda = 1
    muhat <- median(x)
  } else if (islambdatwo) { # lambda = 2
    muhat <- mean(x)
  } else {
    f.opt <- function(mu, lambda, x) {
      sum(abs(x - mu) ^ (lambda - 1.0) * sign(x - mu))
    }
    muhat <- uniroot(f = f.opt, interval = c(min(x), max(x)),
                     lambda = lambda, x = x, 
                     tol = 10 ^ -12, #.Machine$double.eps ^ 0.25,
                     maxiter = 1000, trace = 0)$root
  }
  sighat <- (mean(abs(x - muhat) ^ lambda)) ^ (1.0 / lambda)
  pAPD <- function(q, theta1, theta2, lambda, mu = 0, sigma = 1){
    z <- (q - mu) / sigma
    delta <- 2 * theta1 ^ theta2 * (1 - theta1) ^ theta2 / (theta1 ^ theta2 + (1 - theta1) ^ theta2)
    F <- theta1 * (1 - pgamma(delta / lambda * (-pmin(0, z) / theta1) ^ theta2, shape = 1 / theta2)) +
         (1 - theta1) * pgamma(delta / lambda * (pmax(0, z) / (1 - theta1)) ^ theta2, shape = 1 / theta2)
    return(F) 
  }
  uhat <- pAPD(x, theta1, theta2, lambda, muhat, sighat)
  if (version == 1) { # Z_K
    res <- max(((1:n) - 0.5) * log(((1:n) - 0.5) / (n * uhat)) + 
                 (n - (1:n) + 0.5) * log((n - (1:n) + 0.5) / (n * (1 - uhat))))
  }
  if (version == 2) { # Z_A
    res <- -sum(log(uhat) / (n - (1:n) + 0.5) + log(1 - uhat) / ((1:n) - 0.5))
  }
  if (version == 3) { # Z_C
    res <- sum((log((1 / uhat - 1) / ((n - 0.5) / ((1:n) - 0.75) - 1))) ^ 2)
  }
  return(res)
}

stat102 <- function(x, m = NULL) { # AJ
  n <- length(x)
  if (is.null(m) || (m == 0)) {
    if ((1 <= n) && (n <= 8)) m <- 2
    if ((9 <= n) && (n <= 15)) m <- 3
    if ((16 <= n) && (n <= 25)) m <- 5
    if ((26 <= n) && (n <= 35)) m <- 6
    if ((36 <= n) && (n <= 45)) m <- 7
    if ((46 <= n) && (n <= 60)) m <- 8
    if ((61 <= n) && (n <= 90)) m <- 9
    if ((91 <= n) && (n <= 120)) m <- 10
    if (n >= 121) m <- round((log(n)) ^ 1.5)
  }
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  F <- function(x, mu, sigma) {
    rez <- 0.5 * exp(-abs(x - mu) / sigma)
    rez[x >= mu] <- 1 - rez[x >= mu]
    return(rez)
  }
  res <- 0
  for (j in 1:n) {
    idx1 <- if ((j + m) > n) n else j + m
    idx2 <- if ((j - m) < 1) 1 else j - m
    res <- res + 1 / (F(x[idx1], muhat, sighat) - F(x[idx2], muhat, sighat))
  }
  res <- 2 * m * res / n ^ 2
  return(res)
}

stat103 <- function(x, m0, version, criterion, orthofam) { # CH_{criterion, orthogonal family}^{(\xi)}
  # version = Z,T; criterion = AIC, BIC; orthofam = poly, cos
  n <- length(x)
  x <- sort(x)
  if (length(unique(x)) == 1) stop("All values of x are the same in stat103")
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  sighat <- mean(abs(x - muhat))
  if (abs(sighat) < 10^(-30)) stop("sighat too small in stat103")
  zhat <- (x - muhat) / sighat
  psi <- function(orthofam, j, x) {
    if (orthofam == 1) {
      if (j == 0) rez <- 1.0
      if (j == 1) rez <- 0.707107 * x    
      if (j == 2) rez <- 0.223607 * (-2.0 + x ^ 2.0)
      if (j == 3) rez <- 0.0481125 * (-12.0 * x + x ^ 3.0)
      if (j == 4) rez <- 0.00763274 * (-24.0 + x ^ 4.0 - 33.6 * (-2.0 + x ^ 2.0))
      if (j == 5) rez <- 0.000977577 * (-360.0 * x + x ^ 5.0 - 73.3333 * (-12.0 * x + x ^ 3.0))
      if (j == 6) rez <- 0.000103512 * (-720.0 + x ^ 6.0 - 1944.0 * (-2.0 + x ^ 2.0) - 134.295 * (-24.0 + x ^ 4.0 - 33.6 * (-2.0 + x ^ 2.0)))
      if (j == 7) rez <- 9.44636 * 10.0 ^ -6.0 * (-20160.0 * x + x ^ 7.0 - 7280.0 * (-12.0 * x + x ^ 3.0) - 223.486 * (-360.0 * x + x ^ 5.0 - 73.3333 * (-12.0 * x + x ^ 3.0)))
      if (j == 8) rez <- 7.50613 * 10.0 ^ -7.0 * (-40320.0 + x ^ 8.0 - 177408.0 * (-2.0 + x ^ 2.0) - 20904.2 * (-24.0 + x ^ 4.0 - 33.6 * (-2.0 + x ^ 2.0)) - 343.561 * (-720.0 + x ^ 6.0 - 1944.0 * (-2.0 + x ^ 2.0) - 134.295 * (-24.0 + x ^ 4.0 - 33.6 * (-2.0 + x ^ 2.0))))
      if (j == 9) rez <- 5.32183 * 10.0 ^ -8.0 * (-1.8144 * 10.0 ^ 6.0 * x + x ^ 9.0 - 1.008 * 10.0 ^ 6.0 * (-12.0 * x + x ^ 3.0) - 51546.7 * (-360.0 * x + x ^ 5.0 - 73.3333 * (-12.0 * x + x ^ 3.0)) - 501.94 * (-20160.0 * x + x ^ 7.0 - 7280.0 * (-12.0 * x + x ^ 3.0) - 223.486 * (-360.0 * x + x ^ 5.0 - 73.3333 * (-12.0 * x + x ^ 3.0))))
      if (j == 10) rez <- 3.38415 * 10.0 ^ -9.0 * (-3.6288 * 10.0 ^ 6.0 + x ^ 10.0 - 2.35872 * 10.0 ^ 7.0 * (-2.0 + x ^ 2.0) - 4.15039 * 10.0 ^ 6.0 * (-24.0 + x ^ 4.0 - 33.6 * (-2.0 + x ^ 2.0)) - 111817.0 * (-720.0 + x ^ 6.0 - 1944.0 * (-2.0 + x ^ 2.0) - 134.295 * (-24.0 + x ^ 4.0 - 33.6 * (-2.0 + x ^ 2.0))) - 700.875 * (-40320.0 + x ^ 8.0 - 177408.0 * (-2.0 + x ^ 2.0) - 20904.2 * (-24.0 + x ^ 4.0 - 33.6 * (-2.0 + x ^ 2.0)) - 343.561 * (-720.0 + x ^ 6.0 - 1944.0 * (-2.0 + x ^ 2.0) - 134.295 * (-24.0 + x ^ 4.0 - 33.6 * (-2.0 + x ^ 2.0)))))
    }
    if (orthofam == 2) {
      if (j == 0) rez <- 1.0
      if (j == 1) rez <- 1.69031 * (-0.5 + cos(x))
      if (j == 2) rez <- 1.63272 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x))
      if (j == 3) rez <- 1.63361 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x))
      if (j == 4) rez <- 1.63516 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x))
      if (j == 5) rez <- 1.63602 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x))
      if (j == 6) rez <- 1.63652 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x))
      if (j == 7) rez <- 1.63683 * (-0.02 - 0.0320166 * (-0.5 + cos(x)) - 0.0397866 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0579348 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0983732 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.203683 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.53673 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) +  cos(6.0 * x)) + cos(7.0 * x))
      if (j == 8) rez <- 1.63705 * (-0.0153846 - 0.024015 * (-0.5 + cos(x)) - 0.0282147 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0377136 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.056539 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.0974017 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.203037 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x)) - 0.536415 * (-0.02 - 0.0320166 * (-0.5 + cos(x)) - 0.0397866 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0579348 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0983732 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.203683 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.53673 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x)) + cos(7.0 * x)) + cos(8.0 * x))
      if (j == 9) rez <- 1.6372 * (-0.0121951 - 0.0187007 * (-0.5 + cos(x)) - 0.0211108 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0265633 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0365565 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.0556919 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.0967804 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x)) - 0.202607 * (-0.02 - 0.0320166 * (-0.5 + cos(x)) - 0.0397866 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0579348 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0983732 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.203683 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.53673 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x)) + cos(7.0 * x)) - 0.536196 * (-0.0153846 - 0.024015 * (-0.5 + cos(x)) - 0.0282147 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0377136 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.056539 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.0974017 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.203037 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x)) - 0.536415 * (-0.02 - 0.0320166 * (-0.5 + cos(x)) - 0.0397866 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0579348 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0983732 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.203683 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.53673 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) -  0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) -  0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) +  cos(5.0 * x)) + cos(6.0 * x)) + cos(7.0 * x)) + cos(8.0 * x)) + cos(9.0 * x))
      if (j == 10) rez <- 1.63731 * (-0.00990099 - 0.0149869 * (-0.5 + cos(x)) - 0.0164291 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0197798 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0256038 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.0358319 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.0551364 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x)) - 0.0963581 * (-0.02 - 0.0320166 * (-0.5 + cos(x)) - 0.0397866 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0579348 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0983732 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.203683 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.53673 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x)) + cos(7.0 * x)) - 0.202305 * (-0.0153846 - 0.024015 * (-0.5 + cos(x)) - 0.0282147 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0377136 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.056539 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.0974017 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.203037 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x)) - 0.536415 * (-0.02 - 0.0320166 * (-0.5 + cos(x)) - 0.0397866 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0579348 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0983732 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.203683 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.53673 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) -  0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) -  0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) +  cos(5.0 * x)) + cos(6.0 * x)) + cos(7.0 * x)) + cos(8.0 * x)) - 0.536039 * (-0.0121951 - 0.0187007 * (-0.5 + cos(x)) - 0.0211108 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0265633 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0365565 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.0556919 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.0967804 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) + cos(6.0 * x)) - 0.202607 * (-0.02 - 0.0320166 * (-0.5 + cos(x)) - 0.0397866 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0579348 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.0983732 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.203683 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.53673 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) -  0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) +  cos(5.0 * x)) + cos(6.0 * x)) + cos(7.0 * x)) - 0.536196 * (-0.0153846 - 0.024015 * (-0.5 + cos(x)) - 0.0282147 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.0377136 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) - 0.056539 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.0974017 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) - 0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) + cos(5.0 * x)) - 0.203037 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) - 0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) - 0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) -  0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) -  0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) -  0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) +  cos(5.0 * x)) + cos(6.0 * x)) - 0.536415 * (-0.02 - 0.0320166 * (-0.5 + cos(x)) - 0.0397866 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) - 0.0579348 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) + cos(3.0 * x)) - 0.0983732 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) -  0.203683 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) -  0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) -  0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) +  cos(5.0 * x)) - 0.53673 * (-0.027027 - 0.0449064 * (-0.5 + cos(x)) -  0.0605593 * (-0.2 - 0.571429 * (-0.5 + cos(x)) +  cos(2.0 * x)) -  0.100032 * (-0.1 - 0.226891 * (-0.5 + cos(x)) -  0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) -  0.20472 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) -  0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) - 0.537212 * (-0.0384615 - 0.0676986 * (-0.5 + cos(x)) -  0.103347 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) -  0.206563 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) -  0.538006 * (-0.0588235 - 0.113769 * (-0.5 + cos(x)) - 0.210583 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) - 0.539487 * (-0.1 - 0.226891 * (-0.5 + cos(x)) - 0.543424 * (-0.2 - 0.571429 * (-0.5 + cos(x)) + cos(2.0 * x)) + cos(3.0 * x)) + cos(4.0 * x)) +  cos(5.0 * x)) + cos(6.0 * x)) + cos(7.0 * x)) + cos(8.0 * x)) + cos(9.0 * x)) + cos(10.0 * x))
    }
    return(rez)
  }
  psijbar <- function(orthofam, j, z) mean(psi(orthofam, j, z))
  phi <- function(z) 0.5 * exp(-abs(z))
  integrand.fn <- function(z, orthofam, S, a) {
    p <- length(z)
    rez <- rep(NA, p)
    for (i in 1:p) {
      rez[i] <- 0
      for (j in 1:length(S)) {
        rez[i] <- rez[i] + a[j] * psi(orthofam, S[j], z[i])
      }
      rez[i] <- exp(rez[i]) * phi(z[i])
    }
    return(rez)
  }
  integrandprime.fn <- function(z, orthofam, S, a, j) {
    p <- length(z)
    rez <- rep(NA, p)
    for (i in 1:p) {
      rez[i] <- 0
      for (k in 1:length(S)) {
        rez[i] <- rez[i] + a[k] * psi(orthofam, S[k], z[i])
      }
      rez[i] <- psi(orthofam, S[j], z[i]) * exp(rez[i]) * phi(z[i])
    }
    return(rez)
  }
  cs <- function(a, orthofam, S) {
    tmp <- try(integrate(integrand.fn, lower = -Inf, upper = Inf, orthofam = orthofam, S = S, a = a)$value, silent = TRUE)
    if (class(tmp) == "try-error" || abs(tmp) > 10.0 ^ 30) return(.Machine$double.xmax) else return(tmp)
    return(tmp)
  }
  csprime <- function(a, orthofam, S, j) {
    tmp <- try(integrate(integrandprime.fn, lower = -Inf, upper = Inf, orthofam = orthofam, S = S, a = a, j = j)$value, silent = TRUE)
    if (class(tmp) == "try-error" || abs(tmp) > 10.0 ^ 30) return(.Machine$double.xmax) else return(tmp)
    return(tmp)
  }
  ahatS <- function(orthofam, S, z) {
    lower <- rep(-Inf, length(S))
    upper <- rep(Inf, length(S))
    if (orthofam == 1) { # Recall that we excluded the sets S with max(S) odd unless S = {1}
      if (max(S) > 1) {
        epsilon <- 0.001
        upper[which.max(S)] <- 0 - epsilon # a_j < 0
      } else { # i.e., S = {1}
        epsilon2 <- 0.00001
        lower <- -(sqrt(2) - epsilon2)
        upper <- (sqrt(2) - epsilon2)
      }
    }
    return(optim(par = rep(0, length(S)), fn = minusCHZwon, gr = minusgrCHZwon, orthofam = orthofam, S = S, z = z, method = "L-BFGS-B", lower = lower, upper = upper)$par)
  }
  CHZwon <- function(a, orthofam, S, z) {
    rez <- 0
    for (j in 1:length(S)) {
      rez <- rez + a[j] * psijbar(orthofam, S[j], z)
    }
    tmp <- cs(a, orthofam, S)
    if ((abs(tmp) == .Machine$double.xmax) || (abs(rez) == .Machine$double.xmax)) return(-10 ^ 30) 
    rez <- rez - log(tmp)
    return(rez)
  }
  minusCHZwon <- function(a, orthofam, S, z) {
    return(-CHZwon(a, orthofam, S, z))
  }
  minusgrCHZwon <- function(a, orthofam, S, z) {
    tmp <- cs(a, orthofam, S)
    if (tmp == .Machine$double.xmax) return(rep(-10 ^ 30, length(S)))
    rez <- rep(NA, length(S))
    for (j in 1:length(S)) {
      tmp2 <- csprime(a, orthofam, S, j)
      tmp3 = -psijbar(orthofam, S[j], z) + tmp2 / tmp
      if (abs(tmp3) == .Machine$double.xmax) rez[j] <- -10 ^ 30 else rez[j] <- tmp3
    }
    return(rez)
  }
  CHZ <- function(orthofam, S, z)  {
    ahat <- ahatS(orthofam, S, z)
    return(2 * length(z) * CHZwon(ahat, orthofam, S, z))
  }
  CHT <- function(orthofam, S, z)  {
    rez <- 0
    for (j in S) {
      tmp <- psijbar(orthofam, j, z)
      rez <- rez + tmp ^ 2
    }
    rez <- length(z) * rez
    return(rez)
  }

  vec <- 1:m0
#  if (orthofam == 1) vec <- c(1, vec[!(vec %% 2)]) # remove odd degrees
  out  <-  list(vec[c()])
  for (i in seq_along(vec)) {
    out  <-  c(out, lapply(out[lengths(out) < m0], function(y) c(y, (vec)[i])))
  }
  allS <- out[-1] # All non-empty subsets of {1, ..., m0}
  #allS <- out # All subsets of {1, ..., m0}

#  Sstar <- allS[[1]]
#  chbest <- if (version == 1) CHZ(orthofam, Sstar, zhat) else if (version == 2) CHT(orthofam, Sstar, zhat)
#  if (criterion == 1) {
#    AICmax <- chbest - 2 * length(Sstar)
#  } else if (criterion == 2) {
#    BICmax <- chbest - log(n) * length(Sstar)
#  }

  AICmax <- BICmax <- 0 # empty set
  chbest <- 0
  Sstar <- "empty"
  
  #if (length(allS) > 1) {
    for (j in 1:length(allS)) {
      S <- allS[[j]]
      if (orthofam == 1) {
        if (max(S) > 1) {
          if (!(max(S) %% 2 == 0)) next # Do not keep S when max(S) is odd
        }
      }
      ch <- if (version == 1) CHZ(orthofam, S, zhat) else if (version == 2) CHT(orthofam, S, zhat)
      if (criterion == 1) {
        AIC <- ch - 2 * length(S)
        if (AIC > AICmax) {AICmax <- AIC; Sstar <- S; chbest <- ch}
      } else if (criterion == 2) {
        BIC <- ch - log(n) * length(S)
        if (BIC > BICmax) {BICmax <- BIC; Sstar <- S; chbest <- ch}
      } 
    }
  #}
#  cat("SStar: ") ; cat(Sstar); cat(" ")
  res <- chbest
  return(res)
}

stat104 <- function(x) { # BRT_{3,4}
  return(6 * stat36(x) ^ 2 / 7 + 149 * stat37(x) ^ 2 / 165)
}

stat105 <- function(x) { # RB_{3,4}
  return(stat55(x) ^ 2 + stat56(x) ^ 2)
}

stat106 <- function(x, lambda = 1) { # DLO(lambda)
  n <- length(x)
  
  if (lambda < (1.0 - 10 ^ -15)) stop("lambda should be >=1 in stat106")
  invlambda <- 1.0 / lambda
  beta <- 1.0 + invlambda
  gamma <- -digamma(1.0)
  onemgamma <- 1.0 - gamma
  islambdaone <- abs(lambda - 1.0) < 10 ^ -15
  islambdatwo <- abs(lambda - 2.0) < 10 ^ -15
  islambdaonephalf <- abs(lambda - 1.5) < 10 ^ -15
  islambdatwophalf <- abs(lambda - 2.5) < 10 ^ -15
  islambdathree <- abs(lambda - 3.0) < 10 ^ -15
  onesixteenth <- 1.0 / 16.0
  ctepi <- pi ^ 2.0 / 3.0 - 3.0
  ctepibis <- 3.0 - 8.0 / pi
  ctepibisbis <- (3.0 * pi ^ 2.0 - 28.0) / 8.0
  ctegamma <- 0.5 * (2.0 - log(2.0) - gamma)
  ctebeta <- invlambda * beta * trigamma(beta) - invlambda
  ctelambda <- invlambda * (lambda + log(lambda) + digamma(invlambda))
  ctelambdabis <- 1.0 + lambda - lambda ^ 2.0 / (gamma(2.0 - invlambda) * gamma(invlambda))  
  # Equ. (3.14)
  f106 <- function(mu, lambda, x) {
    sum(abs(x - mu) ^ (lambda - 1.0) * sign(x - mu))
  }
  if (islambdaone) { # lambda = 1
    muhat <- median(x)
  } else if (islambdatwo) { # lambda = 2
    muhat <- mean(x)
  } else {
    muhat <- uniroot(f = f106, interval = c(min(x), max(x)),
                     lambda = lambda, x = x, 
                     tol = 10 ^ -12, #.Machine$double.eps ^ 0.25,
                     maxiter = 1000, trace = 0)$root
  }
  sighat <- (mean(abs(x - muhat) ^ lambda)) ^ invlambda

  # Equ. (4.14)
  y <- (x - muhat) / sighat
  absy <- abs(y)
  plafond <- 10 ^ (300 / lambda)
  absy[absy > plafond] <- plafond
  Blambda <- mean((absy) ^ lambda * sign(y))
  Klambda <- sum((absy[absy > 10^-16]) ^ lambda * log(absy[absy > 10^-16])) / n

  if (islambdaone) { # lambda = 1
    # Equ. (4.16)
    denom1.x <- 1.0
    num2.x <- onemgamma
    denom2.x <- ctepi
  } else if (islambdatwo) { # lambda = 2
    # Equ. (4.17)
    denom1.x <- ctepibis
    num2.x <- 0.5 * (2.0 - log(2.0) - gamma)
    denom2.x <- ctepibisbis
  } else {
    denom1.x <- ctelambdabis
    denom2.x <- invlambda * (beta * trigamma(beta) - 1.0)
    num2.x <- ctelambda
  }
  # Equ. (4.15)
  term1.x <- n * Blambda ^ 2.0 / denom1.x
  term2.x <- n * (Klambda - num2.x) ^ 2.0 / denom2.x
  XlambdaAPDstar <- term1.x + term2.x # X_{\lambda}^{APD*}


  if (islambdaone) { # lambda = 1
    if (n %% 2 == 0) { # n even
      alpha1 <- 1.06
      c1 <- -1.856
      alpha2 <- 1.01
      c2 <- -0.422
      alpha3 <- 0.92
      c3 <- -1.950
      alpha4 <- 2.3
      c4 <- 39.349
    } else { # n odd
      alpha1 <- 1.03
      c1 <- -0.281
      alpha2 <- 0.86
      c2 <- -0.198
      alpha3 <- 1.04
      c3 <- -3.827
      alpha4 <- 1.0
      c4 <- 0
    }
  } else if (islambdatwo) { # lambda = 2
    alpha1 <- 0.99
    c1 <- -1.890
    alpha2 <- 1.00
    c2 <- -0.788
    alpha3 <- 1.05
    c3 <- -9.327
    alpha4 <- 1.4
    c4 <- 14.208
  } else if (islambdaonephalf) { # lambda = 1.5
    alpha1 <- 0.99
    c1 <- -0.952
    alpha2 <- 0.99
    c2 <- -0.637
    alpha3 <- 0.55
    c3 <- -3.488
    alpha4 <- 0.5
    c4 <- 2.434
  } else if (islambdatwophalf) { # lambda = 2.5
    alpha1 <- 0.99
    c1 <- -2.981
    alpha2 <- 0.99
    c2 <- -0.844
    alpha3 <- 1.10
    c3 <- -23.104
    alpha4 <- 1.3
    c4 <- 30.028
  } else if (islambdathree) { # lambda = 3
    alpha1 <- 0.97
    c1 <- -3.855
    alpha2 <- 0.98
    c2 <- -0.880
    alpha3 <- 1.14
    c3 <- -95.743
    alpha4 <- 1.2
    c4 <- 103.871
  }
  

  # Equ. (5.1)
  NKlambda <- max(0, Klambda - 0.5 * lambda * Blambda ^ 2)
  if (islambdaone) { # lambda = 1
    if (n %% 2 == 0) { # n even
      # Equ. (5.9) and (5.10)
      VarBlambda <- (1.0 - 1.856 / n ^ 1.06) / n
      VarNKlambdaonefourth <- (onesixteenth * onemgamma ^ -1.5 * ctepi * (1.0 - 1.950 / n ^ 0.92 + 39.349 / n ^ 2.3)) / n
      ENKlambdaonefourth <- onemgamma ^ 0.25 * (1.0 - 0.422 / n ^ 1.01)
    } else { # n odd
      # Equ. (5.11) and (5.12)
      VarBlambda <- (1.0 - 0.281 / n ^ 1.03) / n
      VarNKlambdaonefourth <- (onesixteenth * onemgamma ^ -1.5 * ctepi * (1.0 - 3.827 / n ^ 1.04)) / n
      ENKlambdaonefourth <- onemgamma ^ 0.25 * (1.0 - 0.198 / n ^ 0.86)
    }
  } else if (islambdatwo) { # lambda = 2
    # Equ. (5.13) and (5.14)
    VarBlambda <- ctepibis * (1.0 - 1.890 / n ^ 0.99) / n
    VarNKlambdaonefourth <- (onesixteenth * ctegamma ^ -1.5 * ctepibisbis * (1.0 - 9.327 / n ^ 1.05 + 14.208 / n ^ 1.4)) / n
    ENKlambdaonefourth <- ctegamma ^ 0.25 * (1.0 - 0.788 / n)
  } else {
    # Equ. (5.4) and (5.5) and (5.6)
    VarBlambda <- ctelambdabis * (1.0 + c1 / n ^ alpha1) / n
    tmplambda <- ctelambda
    ENKlambdaonefourth <- tmplambda ^ 0.25 * (1.0 + c2 / n ^ alpha2)
    VarNKlambdaonefourth <- onesixteenth * tmplambda ^ (-1.5) * invlambda * (beta * trigamma(beta) - 1.0) * (1.0 + c3 / n ^ alpha3 + c4 / n ^ alpha4) / n
  }
  

  if (any(islambdaone, islambdatwo, islambdaonephalf, islambdatwophalf, islambdathree)) {
    # Equ. (5.7) and (5.8) and (5.15)
    ZBlambda <- Blambda / sqrt(VarBlambda)
    ZNKlambdaonefourth <- (NKlambda ^ 0.25 - ENKlambdaonefourth) / sqrt(VarNKlambdaonefourth)
    XlambdaAPD <- ZBlambda ^ 2 + ZNKlambdaonefourth ^ 2
  } else {
    # These are not defined because the c_j's and alpha_j's were not computed
    ZBlambda <- NA
    ZNKlambdaonefourth <- NA
    XlambdaAPD <- NA
  }

  if (islambdaone) { # lambda = 1
    # Equ. (6.3)
    ZlambdaEPDstar <- sqrt(n) * (Klambda - onemgamma) / sqrt(ctepi)
  } else if (islambdatwo) { #lambda = 2
    # Equ. (6.4)
    ZlambdaEPDstar <- sqrt(n) * (Klambda - ctegamma) / sqrt(ctepibisbis)
  } else {
    # Equ. (6.2)
    ZlambdaEPDstar <- sqrt(n) * (Klambda - ctelambda) / sqrt(ctebeta)
  }
  # Equ. (6.19)
  mutilde <- min(x)
  sigmatilde <- (mean((x - mutilde) ^ lambda)) ^ invlambda
  ytilde <- (x - mutilde) / sigmatilde
  Klambdatilde <- sum((abs(ytilde)[abs(ytilde) > 10^-16]) ^ lambda * log(abs(ytilde)[abs(ytilde) > 10^-16])) / n
  ZlambdaHEPDstar <- sqrt(n) * (Klambdatilde - ctelambda) / sqrt(ctebeta)
   
  return(list(XlambdaAPD = XlambdaAPD, pvalXlambdaAPD = pchisq(XlambdaAPD, df = 2, lower.tail = FALSE),
              ZNKlambdaonefourth = ZNKlambdaonefourth, pvalZNKlambdaonefourth = 2 * pnorm(abs(ZNKlambdaonefourth), lower.tail = FALSE),
              XlambdaAPDstar = XlambdaAPDstar, pvalXlambdaAPDstar = pchisq(XlambdaAPDstar, df = 2, lower.tail = FALSE), 
              ZBlambda = ZBlambda, pvalZBlambda = 2 * pnorm(abs(ZBlambda), lower.tail = FALSE),
              ZlambdaEPDstar = ZlambdaEPDstar, pvalZlambdaEPDstar = 2 * pnorm(abs(ZlambdaEPDstar), lower.tail = FALSE),
              ZlambdaHEPDstar = ZlambdaHEPDstar, pvalZlambdaHEPDstar = 2 * pnorm(abs(ZlambdaHEPDstar), lower.tail = FALSE)))
}




stat107 <- function(x) { # KP for Kozubowski / Panorska test
  n <- length(x)
  x <- sort(x)
  muhat <- if ((n %% 2) == 0) 0.5 * (x[n / 2] + x[n / 2 + 1]) else x[(n + 1) / 2]
  tmp2 <- sum(pmax(x - muhat, 0))
  if (abs(tmp2) < 0.0000000001) return(1.0) else {
    Khat <- sum(pmax(-(x - muhat), 0)) / tmp2
    return (n * (2 - (1 + sqrt(Khat)) ^ 2 / (1 + Khat)))
  }
}



stat108 <- function(x) { # SD for Subramanian / Dixit test
  n <- length(x)
  x <- sort(x)
  .hsm <- function(x) {
    while (length(x) >= 4) {
      nx <- length(x)
      k <- ceiling(0.5 * nx) - 1
      inf <- x[1:(nx - k)]
      sup <- x[(k + 1):nx]
      diffs <- sup - inf
      i <- which(diffs == min(diffs))
      if(length(i) > 1) i <- mean(i) 
      if (diffs[i] == 0) {
        x <- x[i]
      } else {
        x <- x[i:(i + k)]
      }
    }
    if (length(x) == 3) {
      z <- 2 * x[2] - x[1] - x[3]
      M <- switch(as.character(sign(z)), 
                  "-1" =  mean(x[1:2]), 
                  "1" = mean(x[2:3]), 
                  "0" = x[2])
    } else {
      M <- mean(x)
    }
    M
  }
  theta <- .hsm(x) # nonparametric estimation of the mode (see the source code below)
  n1 <- sum(x <= theta)
  u <- sum (x[n1] - x[1:n1])
  v <- ifelse(n1 < n, sum (x[(n1 + 1):n] - x[n1 + 1]), 0)
  SD <- ifelse(u > 0 | v > 0, u / (u + v), 0.5)
  return (SD)
}

