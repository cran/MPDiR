\name{ConfRates}
\alias{ConfRates}
\docType{data}
\title{ Confidence Ratings for Detection of Luminance Patch}
\description{
  Data from a detection task in which the observer gave ratings (1--5) indicating his confidence that a signal was presented on that trial.  The signal was a 3 deg square patch of luminance, fixated centrally, embedded in a 10 deg field of random luminance pixels.  The signal was 64 x 64 pixels.  During the experiment, 6 luminance contrast levels were randomly presented and at an equal rate, randomly interspersed catch-trials in which the signal was absent.
}
\usage{data(ConfRates)}
\format{
  A data frame with 196 observations on the following 3 variables.
  \describe{
    \item{\code{Contrast}}{a numeric vector indicating the luminance contrast of the signal, 0 is a catch trial.}
    \item{\code{Rating}}{a numeric vector indicating the confidence rating of the observer as to the presence of the signal}
    \item{\code{Stim}}{a logical vector indicating whether the signal was present or not on the trial}
  }
}

\source{
 Data kindly provided by Fr\'ed\'eric Devinck.
}
\references{
 Devinck, F. (2003) 
\emph{Les traitements visuels chez l'homme : strat\'egies de classification de la forme}.
Universit\'e Lyon 2, Doctoral Thesis.
}
\examples{
data(ConfRates)
}
\keyword{datasets}
