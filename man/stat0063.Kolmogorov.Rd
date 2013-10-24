\name{stat0063.Kolmogorov}
\alias{stat0063.Kolmogorov}
\alias{stat0063}

\title{The Kolmogorov test for uniformity}

\description{The Kolmogorov test \code{D_n} for uniformity is used

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

Kolmogorov, A. N. (1933), Sulla determinazione empirica di una legge di distibuziane, 
\emph{Giornale dell'Istituta Italiano degli Attuari}, \bold{4}, 83--91.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Uniformity.tests}} for other goodness-of-fit tests for uniformity.
}
\keyword{statistic, test, uniformity, Kolmogorov}
