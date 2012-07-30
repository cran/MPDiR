\name{thresh.est}
\alias{thresh.est}
\title{ Estimate Threshold Level from Model Object }
\description{
 Estimates stimulus level at which probability of a response is some fixed value given a model object for a probit model.
}
\usage{
thresh.est(p, obj)
}

\arguments{
  \item{p}{numeric vector giving the performance level for the desired threshold }
  \item{obj}{ object inheriting from class \sQuote{glm} from a fit with a binomial family. }
}
\details{
  Calculates the stimulus level(s) for a probit model of a psychometric function for given levels of performance.  For a more general version, see \link[MASS]{dose.p} in pkg{MASS}.
}
\value{
 Vector of numeric indicating the stimulus levels for the given performance levels.
 }
\author{Kenneth Knoblauch}
\seealso{See also \code{\link[MASS]{dose.p}} }
\examples{
data(HSP)
SHR1.glm <- glm(p/100 ~ log(Q), binomial(probit), HSP,
	subset = Obs == "SH" & Run == "R1")
thresh.est(1:3/4, SHR1.glm)
}

\keyword{univar}
\keyword{models}
