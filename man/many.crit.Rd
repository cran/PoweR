\name{many.crit}
\alias{many.crit}
\title{Computation of critical values for several test statistics}
\description{Computation of critical values for several test statistics,
  several n values, and several level values, for a given distribution}
\usage{
many.crit(law.index,stat.indices,M = 10^3,vectn = c(20,50,100),levels = c(0.05,0.1),
          alter = create.alter(stat.indices),law.pars = NULL,parstats = NULL,model = NULL)
}
\arguments{
  \item{law.index }{law index as given by function \code{\link{getindex}}. length(\code{law.index})=1}
  \item{stat.indices }{vector of statistic indices as given by function \code{\link{getindex}}.}
  \item{M }{number of Monte Carlo repetitions to use.}
  \item{vectn }{vector of number of observations for the samples to be generated.}
  \item{levels }{vector of required level values.}
  \item{alter }{named-list with type of test for each statistical test:
    alter[["statj"]]=0, 1 ,2, 3 or 4; for each \eqn{j} in \code{stat.indices} (0: two.sided=bilateral, 1: less=unilateral, 2:
  greater=unilateral, 3: bilateral test that rejects H0 only for large
  values of the test statistic, 4: bilateral test that rejects H0 only
  for small values of the test statistic)} 
  \item{law.pars }{\code{NULL} or a vector of length at most 4 containing 4 possible
  parameters to generate random values from distribution law(\code{law.pars}[\eqn{j}],\eqn{j<=4})}
  \item{parstats }{named-list of parameter values for each statistic to simulate. 
    The names of the list should be \code{stat}\eqn{j}, \eqn{j} taken in \code{stat.indices}.
    If \code{stat}\eqn{j}=\code{NA}, the default parameter values for
    the test statistic \code{stat}\eqn{j} will be used.}
  \item{model }{NOT IMPLEMENTED YET. If \code{NULL}, no model is used.
  If an integer \eqn{i>0}, the model coded in the C function \code{modele}\eqn{i} is used.
  Else this should be an R function that
    takes three arguments: \code{eps} (vector of \eqn{\epsilon} values), \code{thetavec}
    (vector of \eqn{\theta} values) and \code{xvec} (vector or matrix of \eqn{x} values).This function should take a
    vector of errors, generate observations from a model (with
    parameters \code{thetavec} and values \code{xvec}) based on these errors, then
    compute and return the residuals from the model. See function
    modele1.R in directory inst/doc/ for an example in multiple linear regression.}
}
\value{An object of class \code{critvalues}, which is a list where each element of the list contains a matrix 
  for the corresponding statistic. This column matrices are: \eqn{n} values,
  level values, parameters of the test statistic (\code{NA} if none), left critical values and right critical values).
}
\seealso{See \code{\link{print.critvalues}} for a LaTeX output of the
  results of this function.}
\references{
Lafaye de Micheaux, P. and Tran, V. A. (2014), PoweR: Reproducible
Research Tool to ease Monte-Carlo Power Simulation Studies for
Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{}, --. 
}
\author{P. Lafaye de Micheaux, V. A. Tran}

\examples{
critval <- many.crit(law.index=2,stat.indices=c(10,15),M=10^3,vectn=c(20,50,100),
                     level=c(0.05,0.1),alter=list(stat10=3,stat15=3),law.pars=NULL,
                     parstats=NULL)
print(critval,digits=3,latex.output=FALSE)
}
\keyword{univar}
