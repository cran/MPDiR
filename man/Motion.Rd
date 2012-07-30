\name{Motion}
\alias{Motion}
\docType{data}
\title{ Threshold for First- and Second-order Movement with Age}
\description{
  These data were reported by Thibault et al (2007) (their Figure 4) and provide developmental data on the contrast modulation threshold for detection of a moving stimulus for two kinds of motion, called First-order (or luminance based) and Second-order (contrast based).
}
\usage{data(Motion)}
\format{
  A data frame with 112 observations on the following 5 variables.
  \describe{
    \item{\code{Subject}}{a factor with levels \code{S01},\code{S02}, ..., \code{S70} coding the id of each observer}
    \item{\code{LnAge}}{a numeric vector giving the natural log of age in months}
    \item{\code{Mtype}}{a factor with levels \code{FO} \code{SO} coding the type of stimulus, first-order or second-order}
    \item{\code{Sex}}{a factor with levels \code{f} \code{m} indicating the sex of the observer}
    \item{\code{LnThresh}}{a numeric vector giving the natural log of the percent contrast modulation threshold}
  }
}
\details{
  The data were collected using preferential-looking techniques. Percent contrast varies between 0 and 100. 
  }
\source{
Data kindly provided by F. Vital-Durand.
}
\references{
  Thibault, D., Brosseau-Lachaine, O., Faubert, J. and Vital-Durand, F. (2007) Maturation of the sensitivity for luminance and contrast modulated patterns during development of normal and pathological human children. \emph{Vision Research}, 47, 1561--1569

}
\examples{
data(Motion)
}
\keyword{datasets}
