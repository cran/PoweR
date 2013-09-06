\name{law0007.Uniform}
\alias{law0007.Uniform}
\alias{law0007}

\title{The Uniform Distribution}

\description{Random generation for the Uniform distribution with parameters \code{min} and \code{max}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{min} or \code{max} are not specified they assume the default values of 0 and 1, respectively.

The Uniform distribution has density:
\deqn{ \frac{1}{max - min} }{ f(x) = 1/(max - min) }
for \eqn{min \le x \le max}.
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See function \code{\link{runif}} from package stats. \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(7,10000,law.pars=c(2,9))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Uniform}
