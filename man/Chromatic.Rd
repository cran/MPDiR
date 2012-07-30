\name{Chromatic}
\alias{Chromatic}
\docType{data}
\title{ Chromatic Thresholds across the Life Span }
\description{
 Thresholds for equiluminant chromatic differences as a function of age along 3 axes in the CIE xy chromaticiy diagram. Data were obtained using a spatial, two-alternative forced-choice paradigm.  Responses were obtained from preverbal observers using forced-choice preferential looking.
}
\usage{data(Chromatic)}
\format{
  A data frame with 511 observations on the following 4 variables.
  \describe{
    \item{\code{Log2Age}}{a numeric vector giving the log base 2 age of each observer}
    \item{\code{Age}}{a numeric vector giving the age of each observer in years}
    \item{\code{Thresh}}{a numeric vector giving the threshold modulation in distance along the axis from the \sQuote{white} point in CIE xy space along the axis being tested}
    \item{\code{Axis}}{a factor with levels \code{Deutan} \code{Protan} \code{Tritan} giving the nominal axis in CIE xy space along which the threshold was measured.  The terms refer to the confusion lines through the \sQuote{white} for each of three types of congentially dichromatic observer.  In the CIE xy-plane, these correspond to projections of axes that isolate modulations in each class of human cone photoreceptor.}
  }
}
\details{Equiluminant modulations along the \sQuote{Protan} and \sQuote{Deutan} axes projected in the CIE xy-plane are not cone-isolating.  However, in the experiment, the chromatic modulations were presented mixed with dynamic, random luminance noise that had been shown to mask luminance components of the stimulus, thus favoring detection by chromatic pathways.
 }

\references{
 Knoblauch, K, Vital-Durand, F, Barbur, JL (2001). Variation of chromatic sensitivity across the life span. \emph{Vision Research}., 41, 1:23--36. 
}
\examples{
data(Chromatic)
}
\keyword{datasets}
