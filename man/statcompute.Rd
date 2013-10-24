\name{statcompute}
\alias{statcompute}
\title{Performs a hypothesis test for the given value of statistic.}
\description{Performs the hypothesis test for those added in the package.}
\usage{
statcompute(stat.index, data, levels = c(0.05,0.1), critvalL = NULL,
            critvalR = NULL, alter = 0, stat.pars = NULL)
}
\arguments{
  \item{stat.index }{one statistic index as given by function \code{\link{getindex}}.}
  \item{data }{sample from which to compute the statistic.}
  \item{levels }{vector of desired significance levels for the test.}
  \item{critvalL }{\code{NULL} or vector of left critival values.}
  \item{critvalR }{\code{NULL} or vector of right critival values.}
  \item{alter }{0: two.sided=bilateral, 1: less=unilateral, 2:
  greater=unilateral, 3: bilateral test that rejects H0 only for large
  values of the test statistic, 4: bilateral test that rejects H0 only
  for small values of the test statistic.} 
  \item{stat.pars }{a vector of parameters. 
    If \code{NULL}, the default parameter values for this statistic wil be used.}
}
\value{A list with the test statistic value, the \eqn{p}-value and \code{alter}.
}
\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\examples{
data <- rnorm(50)
# Shapiro-Wilk test:
statcompute(21, data, levels = c(0.05,0.1), critvalL = NULL, critvalR = NULL,
            alter = 0,stat.pars = NULL)
# Identical to:
shapiro.test(data)
}
\keyword{htest}
