\name{polmer}
\alias{polmer}

\title{
Ordinal Regression with Mixed-effects
}
\description{
Regression models are fit responses that are ordered factors with (or without) random effects.

}
\usage{
polmer(formula, data, lnk = "logit", which.lme4 = "lme4.0",  ...)
}
\arguments{
  \item{formula}{a formula object that may include random effect terms using the format of \code{lmer}.
}
  \item{data}{a data frame that includes the terms of the formula object as columns.
}
  \item{lnk}{character specifying a link function, default \sQuote{logit}.)
}
  \item{which.lme4}{character indicating which version of \pkg{lme4} to use.}
  \item{\dots}{additional named arguments passed along to \code{glmer} or \code{glm}.
}
}
\details{This function is a wrapper that calls \code{glmer} from the \pkg{lme4} package if any random effect terms appear in the formula or \code{glm}, if not.  The response term should be of class \sQuote{integer}, as the function will coerce it to \sQuote{ordered}. 
In the formula object, the random effect should be specified as the second level of random effect with
the intercept removed. See the example below.  This is a glitch for the moment.
}
\value{
An object of class \code{mer} or \code{glm} depending on whether or not any random effect terms are included in the formula object.
}
\author{
Kenneth Knoblauch
}


\seealso{
 \code{\link[stats]{glm}}, \code{\link[lme4]{glmer}}, \code{\link[MASS]{polr}}
}
\examples{
	data(Faces)
	if(require(lme4.0, quietly = TRUE)){
	# no random effects specified - calls glm
	Faces.glm <- polmer(SimRating ~ sibs, Faces)
	# random effect of observer - call glmer
	#  GLITCH: must specify random effect as second level of factor (levels of sibs are 0/1)
	
	Faces.glmer <- polmer(SimRating ~ sibs + (sibs1 - 1 | Obs),
		Faces)
	}
}

\keyword{models}
