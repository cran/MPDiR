\name{Faces}
\alias{Faces}
\docType{data}
\title{ Perceived kinship in image pairs of children }
\description{
  Observers rated relatedness of pairs of images of children that were either siblings or not.
}
\usage{data(Faces)}
\format{
  A data frame with 960 observations on the following 6 variables.
  \describe{
    \item{\code{SimRating}}{a numeric vector containing integer values from 0--10 indicating the observers confidence rating that the pair of images were sibilings or not.}
    \item{\code{sibs}}{a factor with levels \code{0} \code{1}, indicating whether the pair of images was, in fact, not related or related as siblings, respectively.}
    \item{\code{agediff}}{a numeric vector indicating the age difference between the face pairs in months.}
    \item{\code{gendiff}}{a factor with levels \code{diff} \code{same} indicating whether the images were of the same of different gender.}
    \item{\code{Obs}}{a factor with levels \code{S1}--\code{S32} coding the identities of individual observers.}
    \item{\code{Image}}{a factor with levels \code{Im1}--\code{Im30} coding the image pairs of the stimulus set.}
  }
}
\details{
Observers used an 11--point scale (0--10) to rate their confidence in whether the image
pairs were of siblings.
}

\references{
 Maloney, L. T., and Dal Martello, M. F. (2006). Kin recognition and the perceived facial similarity of children. \emph{Journal of Vision}, 6(10):4, 1047--1056, \url{http://journalofvision.org/6/10/4/}. 
}
\examples{
data(Faces)
}
\keyword{datasets}
