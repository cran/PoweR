\name{law0008.Student}
\alias{law0008.Student}
\alias{law0008}

\title{The Student t Distribution}

\description{Random generation for the Student t distribution with \code{df} degrees of freedom.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameter.}

\details{If \code{df} is not specified it assumes the default value of 1.

The t distribution with \code{df = k} degrees of freedom has density:
\deqn{ (\sqrt{k\pi})^{-1}\frac{\Gamma\left(\frac{k+1}{2} \right)}{\Gamma\left(\frac{k}{2} \right)}\left(1+\frac{t^2}{k} \right)^{-\frac{k+1}{2}} }{ 
       f(x) = \Gamma((k+1)/2) / (\sqrt(k\pi) \Gamma(k/2)) (1 + x^2/k)^-((k+1)/2) }
for all real \code{x}. It has mean 0 (for \eqn{k > 1}) and variance \eqn{k/(k-2)} (for \eqn{k > 2}). 
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{\code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(8,10000,law.pars=8)
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Student}
