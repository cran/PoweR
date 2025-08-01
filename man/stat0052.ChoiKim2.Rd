\name{stat0052.ChoiKim2}
\alias{stat0052.ChoiKim2}
\alias{stat0052}

\title{The 2nd Choi-Kim test for the Laplace distribution}

\description{The 2nd Choi-Kim test \code{T_{m,n}^{E}} for the Laplace distribution is used

- to compute its statistic and p-value by calling function \code{\link{statcompute}};

- to compute its quantiles by calling function \code{\link{compquant}} or \code{\link{many.crit}};

- to compute its power by calling function \code{\link{powcomp.fast}} or \code{\link{powcomp.easy}}.
}

\details{If \code{m} is not specified it assumes the default value from the Table 4 (Choi and Kim (2006))
which produces the maximum critical values of the test statistic. 
Note that \eqn{m < (n/2)} where \code{n} is the sample size.
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Choi, B. and Kim, K. (2006), Testing goodness-of-fit for Laplace distribution based on maximum entropy, 
\emph{Statistics}, \bold{40}(6), 517--531.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Laplace.tests}} for other goodness-of-fit tests for the Laplace distribution.
}
\keyword{statistic}
\keyword{test}
\keyword{Laplace}
\keyword{Choi-Kim}
