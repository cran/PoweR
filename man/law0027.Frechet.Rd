\name{law0027.Frechet}
\alias{law0027.Frechet}
\alias{law0027}

\title{The Frechet Distribution}

\description{Random generation for the Frechet distribution with parameters \code{mu}, \code{sigma} and \code{alpha}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{mu}, \code{sigma} and \code{alpha} are not specified they assume the default values of 0, 1 and 1, respectively.

The Frechet distribution with parameters \code{mu = }\eqn{\mu}, \code{sigma = }\eqn{\sigma} and \code{alpha = }\eqn{\alpha} has density:
\deqn{ frac{\alpha}{\sigma}\left(\frac{x-\mu}{\sigma}\right)_{+}^{-\alpha-1}\exp\left\{-\left(\frac{x-\mu}{\sigma}\right)^{-\alpha}\right\} }{ 
       f(x) = \alpha/\sigma ((x - \mu)/\sigma)_{+}^{-\alpha - 1} e^-[((x - \mu)/\sigma)^(-\alpha)] }
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(27,10000,law.pars=c(8,6,2))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Frechet}
