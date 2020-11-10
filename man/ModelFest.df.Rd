\name{ModelFest.df}
\alias{ModelFest.df}
\docType{data}
\title{
Contrast Sensitivity of the Observers from the ModelFest Consortium
}
\description{
The data set contains the log contrast sensitivity measurments for the 43 stimuli used in the ModelFest consortium with 4 repeats for each of the 16 observers.
}
\usage{data(ModelFest.df)}
\format{
  A data frame with 2752 observations on the following 3 variables.
  \describe{
    \item{\code{LContSens}}{a numeric vector indicating the measured log contrast sensitivity, inverse of the log contrast threshold}
    \item{\code{Obs}}{a factor with 16 levels, indicating the initials of the observers}
    \item{\code{Stim}}{a factor with 43 levels coding each of the stimuli standardized for this study}
  }
}
\details{
The ModelFest consortium was composed of a group of investigators from 10 laboratories who undertook to obtain measures of contrast sensitivity on a set of standardized achromatic stimuli, varying in spatial configuration, using a standardized set of procedures.
}
\source{
\url{https://jov.arvojournals.org//5/9/6/modelfestbaselinedata.csv}
}
\references{
Carney, T., Tyler, C.W., Watson, A.B., Makous, W., Beutte, B., Chen, C.C., Norcia, A.M., Klein, S.A.(2000) Modelfest: year one results and plans for future years. 
In: B.E. Rogowitz, T.N. Pappas (eds.) \emph{Proceedings of SPIE: Human Vision and Electronic Imaging V}, vol. 3959, pp. 140--151. SPIE 

Watson, A.B., Ahumada, A.J. (2005) A standard model for foveal detection of spatial contrast. \emph{Journal of  Vision} 5, 717--740  

Watson, A. (2006) The spatial standard observer: A human vision model for display inspection. \emph{SID Symposium Digest of Technical Papers}, vol. 37, pp. 1312--1315, SID 
}
\examples{
data(ModelFest.df)
str(ModelFest.df) 
plot(ModelFest.df)
}
\keyword{datasets}
