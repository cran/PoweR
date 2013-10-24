\name{stat0070.Greenwood}
\alias{stat0070.Greenwood}
\alias{stat0070}

\title{The Greenwood test for uniformity}

\description{The Greenwood test \code{G(n)} for uniformity is used

- to compute its statistic and p-value by calling function \code{\link{statcompute}};

- to compute its quantiles by calling function \code{\link{compquant}} or \code{\link{many.crit}};

- to compute its power by calling function \code{\link{powcomp.fast}} or \code{\link{powcomp.easy}}.
}

\details{Note that \code{n} is the sample size.
}
\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software},
\bold{}, --.

Greenwood, M. (1946), The statistical study of infectious diseases, 
\emph{Journal of Royal Statistical Society Series A}, \bold{109}, 85--110.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{\code{\link{Uniformity.tests}} for other goodness-of-fit tests for uniformity.
}
\keyword{statistic, test, uniformity, Greenwood}
