\name{Gabor}
\alias{Gabor}
\docType{data}
\title{ Detection of a 2 cycle/window Gabor temporal modulation}
\description{
  The data record the observations of one observer performing a Yes/No detection task.  The signal was a temporal luminance modulation of a Gaussian spatial profile (\eqn{\sigma}{sigma} = 2.5 deg) on a CRT screen.  The temporal modulation followed a Gabor function with an envelope of \eqn{\sigma}{sigma} = 160 ms and duration of 4\eqn{\sigma}{sigma}.  The carrier was in sine-phase and of frequency 2 cycles/window width.  On each trial random uniform luminance modulations were presented, and on half of the trials, the signal was added in.  The data represent 16 sessions of 224 trials each.
}
\usage{data(Gabor)}
\format{
  A data frame with 114688 observations on the following 4 variables.
  \describe{
    \item{\code{resp}}{a factor with levels \code{H} \code{FA} \code{M} \code{CR}.  These correspond to the response classifications of the observer and correspond to \sQuote{Hit}, \sQuote{False Alarm}, \sQuote{Miss} and \sQuote{Correct Rejection}.}
    \item{\code{time}}{a factor with levels \code{1}, \code{2}, ..., \code{32} indicating the time points of the stimulus coded as factor levels.}
    \item{\code{Time}}{a numeric vector indicating time points of the stimulus in seconds.}
    \item{\code{N}}{a numeric vector of random uniform noise, the random part of the stimulus.  On half of the trials the Gabor signal was added to the noise.}
  }
}
\details{
  More details on the stimulus conditions are given in the reference indicated below.
}
\source{
  J. P. Thomas, and K. Knoblauch (2005) Frequency and phase contributions to the detection of temporal luminance modulation. \emph{Journal of the Optical Society of America A}, 
 22,  2257--2261.}
\examples{
data(Gabor)
str(Gabor) 
}
\keyword{datasets}
