\name{law0036.AsymmetricLaplace}
\alias{law0036.AsymmetricLaplace}
\alias{law0036}

\title{The Asymmetric Laplace Distribution}

\description{Random generation for the Asymmetric Laplace distribution with parameters \code{mu}, \code{b} and \code{k}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}
\details{If \code{mu}, \code{b} or \code{k} are not specified they assume the default values of 0, 1 and 2, respectively.
}
\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See package \code{VGAM}. See \code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(36,10000,law.pars=c(9,2,6))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Asymmetric Laplace}
