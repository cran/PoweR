\name{getnbparstats}
\alias{getnbparstats}
\title{Get numbers of parameters of test statistics.}
\description{Return the default numbers of parameters of the test
  statistics in the package. 
}
\usage{
getnbparstats(stat.indices = NULL)
}
\arguments{
  \item{stat.indices }{if not \code{NULL}, select only the statistics corresponding to
    this vector of indices.}
}
\value{A vector giving the numbers of parameters of test statistics corresponding to
  the vector of indices.
}
\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{
See \code{\link{getnbparlaws}}, \code{\link{getindex}},
\code{\link{law.cstr}}, \code{\link{stat.cstr}}.
}
\examples{
getnbparstats(c(42:53))
}

\keyword{utilities}
