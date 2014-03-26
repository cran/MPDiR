\name{HSP}
\alias{HSP}
\docType{data}
\title{Hecht, Schlaer and Pirenne (1942) psychometric functions}
\description{
  The data frame contains the per cent of correct responses as a function of average number of quanta per flash of light at the cornea and the frequency with which it was seen for 5 psychometric functions reported in Table V of Hecht et al. (1942)
}
\usage{data(HSP)}
\format{
  A data frame with 30 observations on the following 5 variables.
  \describe{
    \item{\code{Q}}{a numeric vector corresponding to the column labelled as \sQuote{No. of quanta} in the original article.}
    \item{\code{p}}{a numeric vector corresponding to the column labelled as \sQuote{Frequency} with a subtitle \sQuote{per cent} in the original article}
    \item{\code{N}}{a numeric vector indicating the number of presentations of the flash}
    \item{\code{Obs}}{a factor with levels \code{MHP} \code{SH} \code{SS} indicating the initials of the observers}
    \item{\code{Run}}{a factor with levels \code{R1} \code{R2} indicating the session or run of the experiment}
  }
}
\details{
  Hecht et al. (1942) reported the data from 5 psychometric functions for 3 observers (the authors), two of which are replications on the same observer.  Each flash intensity was presented 50 times except for the data for observer SH, for which the presentations numbers are 35 and 40, respectively, for the two runs.
}
\source{
 S. Hecht, S. Schlaer and M.H. Pirenne (1942) Energy, quanta and vision. \emph{Journal of General Physiology}, 25, 819--840.
}

\examples{
data(HSP)
lattice::xyplot(p/100 ~ Q | Obs * Run, HSP,  
	type = c("p", "l"), 
	scales = list(x = list(log = TRUE)), 
	as.table = TRUE)
}
\keyword{datasets}
