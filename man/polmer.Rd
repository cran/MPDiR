\name{polmer}
\alias{polmer}

\title{
Ordinal Regression with Mixed-effects
}
\description{
Regression models are fit responses that are ordered factors with (or without) random effects.

}
\usage{
polmer( ...)
}
\arguments{
  \item{\dots}{dummy arguments as the function is deprecated.  See below for where to find code. }
}
\details{This function is a wrapper that calls \code{glmer} from the \pkg{lme4} package if any random effect terms appear in the formula or \code{glm}, if not.  The response term should be of class \sQuote{integer}, as the function will coerce it to \sQuote{ordered}. 
In the formula object, the random effect should be specified as the second level of random effect with
the intercept removed. See the example below.  This is a glitch for the moment.
This function is now deprecated but the code can be found in the scripts directory in the file, \sQuote{Ch9.R}.
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
#	data(Faces)
#	if(require(lme4.0, quietly = TRUE)){
	# no random effects specified - calls glm
#	Faces.glm <- polmer(SimRating ~ sibs, Faces)
	# random effect of observer - call glmer
	#  GLITCH: must specify random effect as second level of factor (levels of sibs are 0/1)
	
#	Faces.glmer <- polmer(SimRating ~ sibs + (sibs1 - 1 | Obs),
#		Faces)
#	}
}

\keyword{models}
