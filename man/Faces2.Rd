\name{Faces2}
\alias{Faces2}
\docType{data}
\title{ Perceived kinship in image pairs of children }
\description{
  Observers classified relatedness of pairs of images of children that were either siblings or not.
}
\usage{data(Faces2)}
\format{
  A data frame with 960 observations on the following 3 variables.
  \describe{
    \item{\code{Resp}}{a factor with levels \code{1} or \code{0} indicating the observers classifications that the pair of images was of siblings or not, respectively.}
    \item{\code{Stim}}{a factor with levels \code{A} \code{P}, indicating whether the pair of images was, in fact, not related or related as siblings, respectively.}
    \item{\code{Obs}}{a factor with levels \code{S1}--\code{S32} coding the identities of individual observers.}
    \item{\code{Image}}{a factor with levels \code{Im1} to \code{Im30} coding the image pair presented on a given trial.}
  }
}
\details{
Observers judged whether pairs of images of infants were siblings or not in a Yes/No
fashion.
}

\references{
  Maloney, L. T., and Dal Martello, M. F. (2006). Kin recognition and the perceived facial similarity of children. \emph{Journal of Vision}, 6(10):4, 1047--1056, \url{http://journalofvision.org/6/10/4/}. 
}
\examples{
data(Faces)
}
\keyword{datasets}
