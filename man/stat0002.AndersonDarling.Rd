\name{stat0002.AndersonDarling}
\alias{stat0002.AndersonDarling}
\alias{stat0002}

\title{The Anderson-Darling test for normality}

\description{The Anderson-Darling test for normality is used

- to compute its statistic and p-value by calling function \code{\link{statcompute}};

- to compute its quantiles by calling function \code{\link{compquant}} or \code{\link{many.crit}};

- to compute its power by calling function \code{\link{powcomp.fast}} or \code{\link{powcomp.easy}}.
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software},
\bold{}, --.

D'Agostino, R.B. and Stephens, M.A. (1986), \emph{Goodness-of-Fit Techniques}, 
Marcel Dekker, New York. (Table 4.9)
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See package \code{nortest}. See \code{\link{Normality.tests}} for other goodness-of-fit tests for normality.
}
\keyword{statistic, test, normality, Anderson-Darling}
