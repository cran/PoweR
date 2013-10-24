########################
# Auxiliary functions: #
########################
indicator <- function(x,a,b) {
  n <- length(x)
  res <- rep(NA,n)
  for (i in 1:n) res[i] <- if(x>b || x<a) 0 else 1
  return(res)
}

######################
# Density functions: #
######################

dlaw1 <- function(x,mu,b) {
# Density of Laplace(mu,b)  
  return(exp(-abs(x-mu)/b)/(2*b))
}

dlaw2 <- function(x,mu,sigma) {
# Density of N(mu,sigma)
  return(exp(-(x-mu)^2/(2*sigma^2))/(sqrt(2*pi)*sigma))
}

dlaw3 <- function(x,l,s) {
# Density of Cauchy(l,s)
  return(1/(pi*s*(1+((x-l)/s)^2)))
}

dlaw4 <- function(x,mu,s) {
# Density of Logistic(mu,s)
  return((1/s)*exp(-(x-mu)/s)*(1+exp(-(x-mu)/s))^(-2))
}

dlaw5 <- function(x,a,b) {
# Density of gamma(a,b)
  return(x^(a-1)*exp(-b*x)/((1/b)^a*gamma(a)))
}

dlaw6 <- function(x,alpha,beta) {
# Density of Beta(alpha,beta)
  return(x^(alpha-1)*(1-x)^(beta-1)*gamma(alpha+beta)/gamma(alpha)/gamma(beta))
}

dlaw7 <- function(x,a,b) {
# Density of Uniform(a,b)
  return(indicator(x,a,b)/(b-a))
}

dlaw8 <- function(x,k) {
# Density of Stuendt(k)
  return((1+x^2/k)^{-(k+1)/2}*gamma((k+1)/2)/gamma(k/2)/sqrt(k*pi))
}

dlaw9 <- function(x,k) {
# Density of Chi-squared(k)
  return(x^(k/2-1)*exp(-x/2)/gamma(k/2)/2^(k/2))
}

dlaw10 <- function(x,mu,sigma) {
# Density of LogNormal(mu,sigma)
  return(exp(-(log(x)-mu)^2/(2*sigma^2))/(x*sigma*sqrt(2*pi)))
}

dlaw11 <- function(x,lambda,k) {
# Density of Weibull(lambda,k)
  return((lambda/k)*(x/k)^(lambda-1)*exp(-(x/k)^lambda))
}

dlaw12 <- function(x,l,b) {
# Density of Shifted Exponential(l,b)
  return(b*exp(-b*(x-l))*indicator(x,l,Inf))
}

dlaw13 <- function(x,j) {
# Density of Power Uniform(j)
  return(x^(-j/(j+1))/(1+j))
}

dlaw14 <- function(x,k,a,b) {
# Density of Average Uniform(k,a,b)

                                        #  n <- length(x)
#  res <- rep(NA,n)
#  for (i in 1:n) {
#    somme =0;
#    for (j in 0:floor(k*(x[i]-a)/(b-a))) {
#      somme = somme + (-1)^j*choose(k,j)*((x[i]-a)/(b-a)-j/k)^(k-1)
#    }
#    res[i] <- somme
#  }
#  return(indicator(x,a,b)*res*k^k/(factorial(k-1)))

  
#  k <- k+1
#  vecj <- 0:k
#  mat <- outer(k*x,vecj,"-") # n x (k+1) 
#  mat <- t(mat)
#  mat[mat<0] <- 0
#  coef <- (-1)^vecj*choose(k,vecj)
#  mat <- mat^(k-1) * coef
#  res <- apply(mat,2,sum)
#  res <- k*res/factorial(k-1)
  # OR ALSO
  n <- length(x)
  res <- rep(NA,n)
  for (i in 1:n) {
    somme = 0
    for (j in 0:k) {
      somme = somme + (-1)^j*choose(k,j)*((x[i]-a)/(b-a)-j/k)^(k-1)*sign((x[i]-a)/(b-a)-j/k)
    }
    res[i] <- 0.5*indicator(x[i],a,b)*somme*k^k/(factorial(k-1))
  }
  return(res)
 
}

dlaw15 <- function(x,j) {
# Density of UUniform(j)
  return((x^(-j/(1+j))+(1-x)^(-j/(1+j)))/(2*(1+j)))
}

dlaw16 <- function(x,j) {
## Density of VUniform(j)
  return( dlaw14(x-0.5,j+1,0,1)*indicator(x,-Inf,1) + dlaw14(x+0.5,j+1,0,1)*indicator(x,0,Inf) )
}

dlaw17 <- function(x,mu,sigma,nu,tau) {
# Density of Johnson SU
  w <- exp(1/tau^2)
  omega <- -nu/tau
  c <- 1/sqrt((w-1)*(w*cosh(2*omega)+1)/2)
  z <- (x-(mu+c*sigma*sqrt(w)*sinh(omega)))/(c*sigma)
  r <- -nu + tau*asinh(z)
  return(exp(-r^2/2)/(c*sigma/tau)/sqrt((z^2+1)*2*pi))
}

#dlaw18 <- function(x,l) {
## Density of Symmetrical Tukey
#  return()
#}

dlaw19 <- function(x,p,m) {
# Density of Location contaminated
  return((p*exp(-(x-m)^2/2)+(1-p)*exp(-x^2/2))/sqrt(2*pi))
}

dlaw20 <- function(x,g,d) {
# Density of Johnson SB
  return(d*exp(-0.5*(g+d*log(x/(1-x)))^2)/(x*(1-x)*sqrt(2*pi)))
}

dlaw21 <- function(x,xi,omega,alpha) {
# Density of Skew Normal
  return((2/omega)*dnorm((x-xi)/omega)*pnorm(alpha*((x-xi)/omega)))
}

dlaw22 <- function(x,p,d) {
# Density of Scale contaminated
  return(((p/d)*exp(-x^2/(2*d^2))+(1-p)*exp(-x^2/2))/sqrt(2*pi))
}

dlaw23 <- function(x,mu,sigma,xi) {
# Density of Generalized Pareto
  res <- (1-xi*(x-mu)/sigma)^((1-xi)/xi)/sigma
  if (xi > 0) {
    res <- indicator(x,0,mu+sigma/xi)*res
  } else if (xi < 0) {
    res <- indicator(x,-mu,Inf)*res
  } else {
    res <- indicator(x,0,Inf) * exp(-(x-mu)/sigma)/sigma
  }
  return(res)
}

dlaw24 <- function(x,mu,sigma,p) {
# Density of Generalized Error Distribution
  return(p*exp(-abs(x-mu)^p/(sigma^p))/(2*gamma(1/p)*sigma))
}

#dlaw25 <- function(x,) {
## Density of Stable
#  return()
#}

dlaw26 <- function(x,mu,sigma) {
# Density of Gumbel
  return(exp(-exp(-(x-mu)/sigma)-(x-mu)/sigma)/sigma)
}

dlaw27 <- function(x,mu,sigma,alpha) {
# Density of Frechet
  res <- alpha*((x-mu)/sigma)^(-alpha-1)*exp(-((x-mu)/sigma)^(-alpha))/sigma
  res[x<=mu] <- 0
  return(res)
}

dlaw28 <- function(x,mu,sigma,xi) {
# Density of
  if (xi == 0) {
    res <- exp(-exp(-(x-mu)/sigma)-(x-mu)/sigma)/sigma
  } else {
    z <- xi*(x-mu)/sigma
    res <- (1+z)^(-1/xi-1)*exp(-(1+z)^(-1/xi))/sigma
    res[z < -1] <- 0
  }
  return(res)
}

dlaw29 <- function(x,alpha) {
# Density of Generalized Arc Sine
  return((sin(pi*alpha)/pi)*x^(-alpha)*(1-x)^(alpha-1))
}

dlaw30 <- function(x,mu,sigma) {
# Density of Folded Normal
  return(indicator(x,0,Inf)*(dnorm(x,mu,sigma)+dnorm(-x,mu,sigma)))
}

dlaw31 <- function(x,p,m,d) {
# Density of Mixture Normal
  return(p*dnorm(x,m,d)+(1-p)*dnorm(x))
}

dlaw32 <- function(x,a,b) {
# Density of Truncated Normal
  return(indicator(x,a,b)*exp(-x^2/2)/sqrt(2*pi)/(pnorm(b)-pnorm(a)))
}

dlaw33 <- function(x,a) {
## Density of Normal with outliers
  return(dnorm(x))
}

dlaw34 <- function(x,th1,th2,th3,crit) {
## Density of Generalized Eponential Power

  n <- length(x)
  z <- y <- res <- rep(NA,n)
  
  calcul_ctenorm_gep <- function(th1,th2,th3) {
    kk <- integrate(dgep,-10,10,1,1,th1,th2,th3)
    return(1/kk$value)
  }
  calcul_sigma_gep <- function(constante,th1,th2,th3) {
    crit <- 0.30
    int <- function(a) {integrate(dgep,0,a,constante,1,th1,th2,th3)$value}
    sigma <- .1
    flag <- 1
    while (flag == 1) {
      f <- int(sigma)
      if (f>crit) flag <- 0 else sigma <- sigma*1.1; 
    }
    return(1/sigma)
  }
  dgep <- function(y,constante,sigma,th1,th2,th3) {
    n <- length(y)
    z <- res <- rep(NA,n)
    for (i in 1:n) {
      z[i] <- abs(y[i]/sigma)
      res[i] <- (constante/sigma)*exp(-0.5*z[i]^th1)*(1+z[i])^(-th2)*(log(exp(1)+z[i]))^(-th3)
    }
    return(res)
  }
  calcul_z0 <- function(constante,sigma,th1,th2,th3,crit=1e-6) {
    int <- function(a) {integrate(dgep,0,a,constante,sigma,th1,th2,th3)$value}
    z0 <- 1
    flag <- 1
    while (flag == 1) {
      f <- int(z0)*2
      if (f>(1-crit)) flag <- 0 else z0 <- z0*1.1;
    }
    return(z0)
  }
  calcul_sup <- function(constante,sigma,z0,th1,th2,th3,critere=1e-6) {
    pas <- z0/100
    range <- seq(0,z0,pas)
    flag <- 1
    while (flag == 1) {
      ff <- dgep(range,constante,sigma,th1,th2,th3)
      pos <- min(which.max(ff))
      sup <- range[pos]
      if (pas<critere) flag <- 0
      range <- seq(max(0,sup-2*pas),sup+2*pas,pas/100)
      pas <- pas/100
    }
    return(dgep(sup,constante,sigma,th1,th2,th3))
  }
  
  constante <- calcul_ctenorm_gep(th1,th2,th3)
  sigma <- calcul_sigma_gep(constante,th1,th2,th3)
  mu <- 0
  for (i in 1:n) {
    y[i] <- x[i]-mu
    z[i] <- abs(y[i]/sigma)
    res[i] <- (constante/sigma)*exp(-.5*z[i]^th1)*(1+z[i])^(-th2)*(log(exp(1)+z[i]))^(-th3)
  }
  return(res)
}

dlaw35 <- function(x,lambda) {
# Density of Exponential
  return(indicator(x,0,Inf)*lambda*exp(-lambda*x))
}

dlaw36 <- function(x,mu,b,k) {
# Density of Asymmetric Laplace
  n <- length(x)
  res <- rep(NA,n)
  for (i in 1:n) {
    if (x[i]<=mu) {
      res[i] <- (sqrt(2)/b)*(k/(1+k^2))*exp(-sqrt(2)*abs(x[i]-mu)/(b*k))
      } else {
        res[i] <- (sqrt(2)/b)*(k/(1+k^2))*exp(-sqrt(2)*k*abs(x[i]-mu)/b)
      }
  }
  return(res)
}

dlaw37 <- function(x,alpha,beta,delta,mu) {
# Density of Normal-inverse Gaussian
  gamma <- sqrt(alpha^2-beta^2)
  return(alpha*delta*besselK(alpha*sqrt(delta^2+(x-mu)^2),1)*exp(delta*gamma+beta*(x-mu))/pi/sqrt(delta^2+(x-mu)^2))
}

dlaw38 <- function(x,theta=0.0,phi=1.0,alpha=0.5,lambda=2.0) {
# Density of Asymmetric Power Distribution
  n <- length(x)
  res <- rep(NA,n)
  delta <- (2*alpha^lambda*(1-alpha)^lambda)/(alpha^lambda+(1-alpha)^lambda)
  for (i in 1:n) {
    if (x[i] <= theta) {
      res[i] <- delta^(1/lambda)*exp(-delta*((abs(x[i]-theta)/phi)^lambda)/(alpha^lambda))/gamma(1+1/lambda)/phi
    } else {
      res[i] <- delta^(1/lambda)*exp(-delta*(((x[i]-theta)/phi)^lambda)/((1-alpha)^lambda))/gamma(1+1/lambda)/phi
    }
  }
  return(res)
}

