\name{panel.psyfun}
\alias{panel.psyfun}

\title{ Panel Function for Adding Psychometric Function Fit to Each Panel }
\description{
Provides a panel function for use with \link[lattice]{xyplot} from the \pkg{lattice} package so that a psychometric function will be fitted and plotted in each panel of a \code{lattice} plot similar to \code{\link[lattice]{panel.lmline}} for fitting a linear regression line.
}
\usage{
panel.psyfun(x, y, n, lnk = "logit", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of covariate values }
  \item{y}{vector of responses, usually proportion of correct responses }
  \item{n}{ vector of integers specifying number of trial for each proportion in \code{y} }
  \item{lnk}{link function to use in fitting \code{glm} with binomial family (default is \dQuote{logit}).}
  \item{\dots}{ Other parameters passed from the formal arguments of the \code{panel} function. }
}

\value{
 Nothing returned.  Function just used for its side-effects of producing a graph within a panel function.}
\author{Kenneth Knoblauch}
\seealso{See also as \code{\link{lattice}}}
\examples{
lattice:::xyplot(Pc ~ Phaseshift | WaveForm + TempFreq + 
	Direction, Vernier, layout = c(4, 2),
	panel = function(x, y, n = 20, ...) {
		lattice:::panel.xyplot(x, y)
		panel.psyfun(x, y, 20, lnk = "probit")
	}
)

}

\keyword{hplot}
\keyword{dplot}
\keyword{aplot}

