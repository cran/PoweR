\name{law0010.LogNormal}
\alias{law0010.LogNormal}
\alias{law0010}

\title{The Log Normal Distribution}

\description{Random generation for the Log Normal distribution whose logarithm 
  has mean equal to \code{meanlog} and standard deviation equal to \code{sdlog}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{meanlog} or \code{sdlog} are not specified they assume the default values of 0 and 1, respectively.

The Log Normal distribution has density:
\deqn{ \frac{1}{x\sigma\sqrt{2\pi}}e^{-\frac{(\ln x-\mu)^2}{2\sigma^2}} }{ 
	   f(x) = 1/(\sqrt(2\pi) \sigma x) e^-((log x - \mu)^2 / (2\sigma^2)) }
where \eqn{\mu} and \eqn{\sigma} are the mean and standard deviation of the logarithm. 
The mean is \eqn{E(X) = exp(\mu + 1/2 \sigma^2)} and the variance is \eqn{Var(X) = exp(2*\mu + \sigma^2)*(exp(\sigma^2) - 1)}.
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{\code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(10,10000,law.pars=c(8,6))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Log Normal}
