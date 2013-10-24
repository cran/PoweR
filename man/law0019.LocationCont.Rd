\name{law0019.LocationCont}
\alias{law0019.LocationCont}
\alias{law0019}

\title{The Location Contaminated Distribution}

\description{Random generation for the Location Contaminated distribution with parameters \code{p} and \code{m}.

This generator is called by function \code{\link{gensample}} to create random variables based on its parameters.}

\details{If \code{p} or \code{m} are not specified they assume the default values of 0.5 and 0, respectively.

The Location Contaminated distribution has density:
\deqn{ \frac{1}{\sqrt{2\pi}}\left[pe^{-\frac{(x-m)^2}{2}}+(1-p)e^{-\frac{x^2}{2}}\right] }{ 
	   f(x) = 1/(\sqrt(2\pi)) (p e^-((x-m)^2/2) + (1-p) e^-(x^2/2)) }
}

\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. _Journal of Statistical Software_, *?*, -. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\seealso{\code{\link{Distributions}} for other standard distributions.
}

\examples{res <- gensample(19,10000,law.pars=c(8,6))
res$law
res$law.pars
mean(res$sample)
sd(res$sample)
}

\keyword{law, Location Contaminated}
