\name{Context}
\alias{Context}
\docType{data}
\title{
Psychometric Functions for the Effects of Context on Contrast Detection}
\description{
Yes-No psychometric functions in which observers reported whether a vertical bar (the target) was present on each trial in the presence of contextual bars.
}
\usage{data(Context)}
\format{
  A data frame with 120 observations on the following 5 variables.
  \describe{
    \item{\code{Obs}}{a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{F} identifying each observer.}
    \item{\code{ContCntr}}{a numeric vector indicating the contrast of the contextual elements.}
    \item{\code{TargCntr}}{a numeric vector indicating the contrast of the target.}
    \item{\code{NumYes}}{an integer vector indicating the number trials on which the observer responded as having seen the target.}
    \item{\code{NumNo}}{an integer vector indicating the number of trials on which the observer responded as not haveing seen the target.}
  }
}
\details{
The number of trials on which the target was detected is reported as a function of its contrast.  The target bar was presented in the context of 6 aligned vertical bars (3 above and 3 below) and the 4 conditions correspond to 4 levels of contrast of the contextual bars, including a 0 contrast condition, i.e., no contextual bars. 
}
\source{
Data kindly provided by Li Zhaoping
}
\references{
Zhaoping L. \& Jingling L. (2008) Filling-in and suppression of visual perception from context---a Bayesian account of perceptual biases by contextual influences, \emph{PLoS Computational Biology} 4(2): e14 doi:10.1371/journal.pcbi.0040014 }
\examples{
data(Context)
Context <- within(Context, Pc <- NumYes/(NumYes + NumNo))
lattice::xyplot(Pc ~ TargCntr | ContCntr, data = Context, 
	groups = Obs,	type = c("l", "p"),
	auto.key = list(space = "right"))
}
\keyword{datasets}
