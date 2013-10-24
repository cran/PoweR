\name{law0009.Chisquared}
\alias{law0009.Chisquared}
\alias{law0009}

\title{The Chi-Squared Distribution}

\description{Random generation for the Chi-squared distribution with \code{df} degrees of freedom.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameter.}

\details{If \code{df} is not specified it assumes the default value of 1.

The Chi-squared distribution with \code{df = k} degrees of freedom has density:
\deqn{ 2^{-k/2}\Gamma(k/2)^{-1}x^{k/2-1}e^{-x/2} }{ 
       f(x) = 1/(2^(k/2) \Gamma(k/2)) x^(k/2-1) e^(-x/2) }
for \eqn{x > 0} and \eqn{k \ge 1}. The mean and variance are \eqn{n} and \eqn{2n}.
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{\code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(9,10000,law.pars=8)
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Chi-Squared}
