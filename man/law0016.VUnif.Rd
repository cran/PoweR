\name{law0016.VUnif}
\alias{law0016.VUnif}
\alias{law0016}

\title{The VUniform Distribution}

\description{Random generation for the VUniform distribution with parameter \code{size}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameter.}

\details{If \code{size} is not specified it assumes the default value of 1.
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 

Quesenberry and Miller (1977), Power studies of some tests for uniformity, \emph{Journal of Statistical Computation and Simulation}, \bold{5}(3), 169--191 (see p. 179)
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{See \code{\link{law0007.Uniform}} for the Uniform
  distribution. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(16,10000,law.pars=8)
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, VUniform}
