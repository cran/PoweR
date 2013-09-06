\name{law0011.Weibull}
\alias{law0011.Weibull}
\alias{law0011}

\title{The Weibull Distribution}

\description{Random generation for the Weibull distribution with parameters \code{shape} and \code{scale}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{shape} or \code{scale} are not specified they assume the default values of 1 and 1, respectively.

The Weibull distribution with \code{shape} parameter \eqn{k} and \code{scale} parameter \eqn{\lambda} has density given by
\deqn{ \frac{k}{\lambda}\left(\frac{x}{\lambda}\right)^{k-1}e^{-(x/\lambda)^k} }{ 
	   f(x) = (k/\lambda) (x/\lambda)^(k-1) e^(-(x/\lambda)^k) }
for \eqn{x > 0}. The cumulative distribution function is \eqn{F(x) = 1 - e^(-(x/\lambda)^k)} on \eqn{x > 0}, 
the mean is \eqn{E(X) = \lambda \Gamma(1 + 1/k)}, and the \eqn{Var(X) = \lambda^2 * (\Gamma(1 + 2/k) - (\Gamma(1 + 1/k))^2)}.
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{\code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(11,10000,law.pars=c(8,6))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Weibull}
