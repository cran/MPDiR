\name{binom.diagnostics}
\alias{binom.diagnostics}
\alias{plot.binom.diag}
\title{
Diagnostics for Binary GLM
}
\description{
Two techniques for evaluating the adequacy of the binary glm model used in \code{mlds}, based on code in Wood (2006).
}
\usage{
binom.diagnostics(obj, nsim = 200, type = "deviance", no.warn = TRUE)

\method{plot}{binom.diag}(x, alpha = 0.025, breaks = "Sturges", ...)
}

\arguments{
  \item{obj}{
Object of class \sQuote{glm} from a fit to data with a binary response.
}
  \item{nsim}{
number of bootstrap simulations to run.
}
  \item{type}{
character indicating type of residuals (defaults to \dQuote{deviance}) to use in the simulations.
}
  \item{no.warn}{
logical indicating whether to suppress warnings from \code{glm}.  Defaults to TRUE.
}
\item{x}{object of class \sQuote{binom.diag}, typically obtained from running \code{binom.diagnostics}}
\item{alpha}{confidence level for envelope of the distribution of residuals}
\item{breaks}{character indicating algorithm for choosing the breaks in the histogram of runs of residuals.  Defaults to \dQuote{Sturges}.  See \code{\link{hist}}.}
\item{...}{additional parameters specifications for the empirical cdf plot}
}
\details{
Wood (2006) describes two diagnostics of the adequacy of a binary glm model based on analyses of residuals (see, p. 115, Exercise 2 and his solution on pp 346-347).  The first one compares the empirical cdf of the deviance residuals to a bootstrapped confidence envelope of the curve.  The second examines the number of runs in the sorted residuals with those expected on the basis of independence in the residuals, again using a resampling based on the models fitted values. The plot method generates two graphs, the first being the empirical cdf and the envelope.  The second is a histogram of the number of runs from the bootstrap procedure with the observed number indicated by a vertical line.  Currently, this only works if the \sQuote{glm} method is used to perform the fit and \emph{not} the \sQuote{optim} method
}
\value{
\code{binom.diagnostics} returns a list of class
\sQuote{binom.diag} with components
  \item{NumRuns }{integer vector giving the number of runs obtained for each simulation}
  \item{resid }{numeric matrix giving the sorted deviance residuals in each column from each simulation}
  \item{Obs.resid }{numeric vector of the sorted observed deviance residuals}
  \item{ObsRuns }{integer giving the observed number of runs in the sorted deviance residuals}
  \item{p }{numeric giving the proportion of runs in the simulation less than the observed value.}

}
\references{
Wood, SN \emph{Generalized Additive Models: An Introduction with R}, Chapman \& Hall/CRC, 2006.
}
\author{
Ken Knoblauch}


\examples{
\dontrun{
# generate psychometric function and binary data
levs <-  10^seq(-2, 0, len = 5)
psyF <-  pnorm(levs, mean = 0.2, sd = 0.2)
Ntrials <- 100
Resp <- rbinom(Ntrials * length(psyF), 1, psyF)
d.df <- data.frame(resp = Resp, levs = levs)
fit.glm <-  glm(resp ~ levs, binomial, d.df)
fit.diag <- binom.diagnostics(fit.glm)
plot(fit.diag)
}
}

\keyword{models}
\keyword{hplot}
