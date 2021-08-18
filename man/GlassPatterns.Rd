\name{GlassPatterns}
\alias{GlassPatterns}
\docType{data}
\title{ Detection of concentric patterns of dipoles in Glass Patterns}
\description{
  The data record the observations of one observer performing a Yes/No detection task.  The signal was a proportion of dot pairs oriented concentrically on virtual circles in the presence of randomly oriented dot pairs.  The proportion of oriented dot pairs is the coherence.  On one-sixth of the trials, the coherence was 0, and from these a false alarm rate can be established.  The color of the dots was varied across trials and is designated by an azimuthal angle in the equiluminant plane.
}
\usage{data(Gabor)}
\format{
  A data frame with 48 observations on the following 5 variables.
  \describe{
    \item{\code{Pyes}}{a numeric indicating the proportion of trials on which the observer responded Yes, the pattern was present.}
    \item{\code{Direction}}{a factor with levels \code{A0}, \code{A45}, ..., \code{A315} indicating the hue directions in the equiluminant plane of the dots.}
    \item{\code{Coherence}}{a numeric vector indicating the proportion of dot pairs falling along virtual concentric circles.}
    \item{\code{N}}{a numeric vector indicating the number of trials upon which the proportion of Yes responses is based.}
    \item{\code{Angle}}{a numeric vector indicating the azimuthal angle in the equiluminant plane of the hue of the dots.}
  }
}
\details{
  The data are unpublished observations from one observer from a larger data set.
}

\examples{
data(GlassPatterns)
str(GlassPatterns) 
}
\keyword{datasets}
