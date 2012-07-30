\name{psyfun.boot}
\alias{psyfun.boot}
\title{Bootstrapping Standard Errors of Psychometric Function Parameters }
\description{
  A function that will run a bootstrap on the estimated parameters of a psychometric function fit given a model object.
}
\usage{
psyfun.boot(obj, N = 100)
}

\arguments{
  \item{obj}{object inheriting from class \sQuote{glm} from a fit of a psychometric function }
  \item{N}{integer indicating number of bootstrap replications. }
}
\details{
   The function computes new binomial responses based on the fitted probabilities of the model object for each bootstrap replication.  A psychometric function is then fit to each one and the fitted coefficients returned as a bootstrap replicate.
}
\value{
 Returns a matrix with one row for each coefficient of the model and one column for each bootstrap replication.
 }
\references{
Maloney, L. T. (1990) Confidence interval for the parameters of psychometric functions.  
\emph{Perception \& Psychophysics}, \bold{47(2)}, 127--134.

Foster, D.H., Bischof, W.F.(1997) Bootstrap estimates of the statistical accuracy of 
thresholds obtained from psychometric functions. \emph{Spatial Vision}, \bold{11(1)}, 135--139.

Treutwein, B., Strasburger, H. (1999) Fitting the psychometric function. \emph{Perception \& 
Psychophysics}, \bold{61(1)}, 87--106.
}

\author{Kenneth Knoblauch}

\examples{
data(HSP)
SHR2 <- subset(HSP, Obs == "SH" & Run == "R2")
SHR2 <- within(SHR2, {
	nyes <- N * p/100
	nno <- N - nyes
	})
SHR2.glm <- glm(cbind(nyes, nno) ~ log(Q), binomial, SHR2)
### For a real problem, set N to 10000 or so
SHR2.boot <- psyfun.boot(SHR2.glm, 10)
}

\keyword{models}
\keyword{nonparametric}
