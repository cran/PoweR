\name{stat.cstr}
\alias{stat.cstr}
\title{
Gives information about a test statistic.
}
\description{
To obtain the name of a test as well as its default number of parameters
and default parameter values.
}
\usage{
stat.cstr(stat.index, stat.pars = NULL)
}
\arguments{
  \item{stat.index }{a single integer value corresponding to the index of a test statistic as
given by function \code{\link{getindex}}.
}
  \item{stat.pars }{vector of the values of the parameters of the test specified in
\code{stat.index}. If \code{NULL}, the default values are used.
}
}
\value{
\item{name }{name of the test.}
\item{nbparams }{default number of parameters of the test.}
\item{law.pars }{values of the parameters}
\item{alter }{0: two.sided=bilateral, 1: less=unilateral, 2:
  greater=unilateral, 3: bilateral test that rejects H0 only for large
  values of the test statistic, 4: bilateral test that rejects H0 only
  for small values of the test statistic.}
}
\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{law.cstr}}, \code{\link{getindex}},
\code{\link{getnbparlaws}}, \code{\link{getnbparstats}}.}
\examples{
stat.cstr(80)
}
\keyword{utilities}
