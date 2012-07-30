\name{lpois}
\alias{lpois}
\alias{lpois1}
\alias{lnorm}

\title{ Likelihoods for Poisson and Gaussian psychometric functions }
\description{
  These functions define likelihoods for binomial models with Gaussian or Poisson psychometric functions.
}
\usage{
lpois(p, d)
lpois1(q, p, d)
lnorm(p, d)
}
\arguments{
  \item{p}{ numeric vector of parameters (in \code{lpois} and \code{lnorm}) or fixed Poisson parameter (\code{lpois2}). }
  \item{d}{ data frame with number of correct (\code{nyes}) and incorrect (\code{nno}) responses and intensity levels (\code{Q}) }
  \item{q}{numeric vector of length 1 giving scale parameter (\code{lpois2}) }
}
\details{
 These functions calculate the log likelihood for binomial models with Gaussian or Poisson psychometric functions and can be used with \link[stats]{optim} to find the best parameters.
}
\value{
 Returns a numeric value giving minus the log likelihood for the given model.}

\author{Kenneth Knoblauch }

\seealso{See \code{\link[stats]{optim}} }
\examples{
data(HSP)
SHR2 <- subset(HSP, Obs == "SH" & Run == "R2")
SHR2 <- within(SHR2, {
	nyes <- N * p/100
	nno <- N - nyes
	})
optim(par = c(5, 0.35), lnorm, d = SHR2)
}
\keyword{models}
