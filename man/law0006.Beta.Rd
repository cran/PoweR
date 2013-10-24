\name{law0006.Beta}
\alias{law0006.Beta}
\alias{law0006}

\title{The Beta Distribution}

\description{Random generation for the Beta distribution with parameters \code{shape1} and \code{shape2}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{shape1} or \code{shape2} are not specified they assume the default values of 1 and 1, respectively.

The Beta distribution with parameters \code{shape1} = \code{a} and \code{shape2} = \code{b} has density:
\deqn{ \frac{\Gamma(\alpha+\beta)}{\Gamma(\alpha)\Gamma(\beta)}x^{\alpha-1}(1-x)^{\beta-1} }{ f(x) = \Gamma(a+b)/(\Gamma(a)\Gamma(b)) x^(a-1) (1-x)^(b-1) }
for \eqn{a > 0, b > 0} and \eqn{0 \le x \le 1} where the boundary values at \eqn{x=0} or \eqn{x=1} are defined as by continuity (as limits). 
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See function \code{\link{rbeta}} from package stats. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(6,10000,law.pars=c(9,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Beta}
