\name{Grue}
\alias{Grue}
\docType{data}
\title{ Distribution of Grue Languages with UV-B Exposure}
\description{
  Data reported by Lindsey and Brown (2002) based on a corpus of 203 languages (only 201 represented in this data set) of the numbers of languages with separate words for \sQuote{blue}
  and \sQuote{green}, a combined word for the two, so-called grue-languages, or the word \sQuote{dark} for the two, as a function of the average annual exposure to UV-B light.
}
\usage{data(Grue)}
\format{
  A data frame with 4 observations on the following 4 variables.
  \describe{
    \item{\code{UV_B}}{a numeric vector indicating average annual UV-B exposure in kiloJoules/m^2}
    \item{\code{blue-green}}{a numeric vector giving the numbers of languages with a separate terms for \sQuote{blue} and \sQuote{green} }
    \item{\code{grue}}{a numeric vector, the numbers of languages with a single term for \sQuote{green/blue}}
    \item{\code{dark}}{a numeric vector, the numbers of languages with a term \sQuote{dark} for green or blue}
  }
}

\source{
Data kindly provided by Angela M. Brown.
 }
\references{
 D. T. Lindsey and A. M. Brown (2002) Color naming and the phototoxic effects of sunlight on the eye. \emph{Psychological Science}, \bold{13}, 506-512.
}
\examples{
data(Grue)
barplot(t(as.matrix(Grue[, -1])), names.arg = round(Grue$UV_B, 1), 
	beside = TRUE, legend.text = c("blue-green", "grue", "dark"),
	col = c("white", "grey", "black"),
	ylim = c(0, 40), cex.lab = 1.5,
	xlab = expression(paste("Average UV-B (kJ/", m^2, ")")),
	ylab = "Numbers of Languages"
	)

}
\keyword{datasets}
