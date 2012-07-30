\name{Vernier}
\alias{Vernier}
\docType{data}
\title{ Psychometric Functions for Vernier Judgement of Phase Difference between Drifting Gratings}
\description{
 A pair of adjacent, horizontal luminance gratings (sine- or square-wave spatial profiles, 0.8 cycles/degree and equal contrast = 0.30) drifting at 2 or 8 Hz.  The gratings were separated by a 30-arcmin lateral gap.  On each trial, the grating pair was drifting either upward or downward and the phase shift between the gratings was chosen at one of 8 values between -50 and 50 degrees.  The observer was instructed to report which grating appeared shifted upward.  The method of constant stimuli was used.
}
\usage{data(Vernier)}
\format{
  A data frame with 64 observations on the following 8 variables.
  \describe{
    \item{\code{Phaseshift}}{a numeric vector giving the Phaseshift between the pair of gratings}
    \item{\code{WaveForm}}{a factor with levels \code{Sine} \code{Square} giving the spatial luminance profile of the gratings}
    \item{\code{TempFreq}}{a factor with levels \code{2} \code{8} indicating the temporal frequency in Hertz of the moving gratings}
    \item{\code{Pc}}{a numeric vector, indicating the proportion of trials on which a particular grating appeared to be shifted upward}
    \item{\code{Direction}}{a factor with levels \code{Downward} \code{Upward} indicating the direction of motion of the grating pair.}
    \item{\code{N}}{a numeric vector indicating the number of trials on which a condition was presented}
    \item{\code{NumUpward}}{a numeric vector indicating the number of trials on which the observer responded that a particular grating was shifted upward}
    \item{\code{NumDownward}}{a numeric vector indicating the number of trials on which the observer reported a particular grating as shifted downward}
  }
}

\source{
 Data kindly provided by H. Sun.
}
\references{
  H. Sun, B. B. Lee and R. C. Baraas Systematic misestimation in a vernier task arising from contrast mismatch, \emph{Visual Neuroscience}, 2008, 25, 365--370.
}
\examples{
data(Vernier)

}
\keyword{datasets}
