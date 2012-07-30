\name{StairCase}
\alias{StairCase}
\docType{data}
\title{ 3-1 Randomly Interleaved Staircase Detection Experiment}
\description{
  Data from one observer performing a 2-alternative forced-choice detection experiment.  A 3-1 staircase procedure was used with two randomly interleaved staircases. 
}
\usage{data(StairCase)}
\format{
  A data frame with 96 observations on the following 4 variables.
  \describe{
    \item{\code{Trial}}{a numeric vector indicating the trial number for the staircase}
    \item{\code{Contrast}}{a numeric vector indicating the contrast of the test stimulus}
    \item{\code{Response}}{a numeric vector taking the values of 1, when the correct interval was detected, and 0 when it was not}
    \item{\code{StairCase}}{a factor with levels \code{S1} \code{S2} indicating to which staircase the trial belonged.}
  }
}
\details{
  The data were obtained using a 3-down, 1-up staircase procedure.  The observer was required to make 3 successive correct judgments before the stimulus contrast was lowered; one incorrect response resulted in an increase in contrast at the next trial.  Initially, the contrast was reduced by a factor of 2 for each sequence of 3 correct responses. After a first incorrect response, successive increments and decrements were by a factor of 1.26.  Each staircase was 48 trials long.  In the actual experiment, the two staircases were run concurrently, with the contrast level chosen by switching from one staircase to the other at random.
  }
\examples{
clr = c("red", "blue")
with(StairCase, plot(Contrast ~ Trial, log = "y",
	col = clr[unclass(StairCase)], pch = 16))
}
\keyword{datasets}
