\name{CorticalCells}
\alias{CorticalCells}
\docType{data}
\title{
Contrast Response Functions from Visual Cortical Cells
}
\description{
The responses in impulses/second of 6 cells recorded in the macaque cortex in visual areas V1 and V2 to a series of contrasts of drifting sinusoidal gratings of optimal location, orientation, size and spatial frequency.
}
\usage{data(CorticalCells)}
\format{
  A data frame with 42 observations on the following 7 variables.
  \describe{
    \item{\code{Contrast}}{numeric vector of the Michelson contrast of the grating}
    \item{\code{Response}}{numeric vector of the mean response of the cell in impulse/second (action potentials)}
    \item{\code{SEM}}{numeric vector of the standard error of the mean of the \code{Response}}
    \item{\code{N}}{numeric vector of the number of responses recorded}
    \item{\code{Type}}{factor with levels \code{CX} \code{SM} corresponding to the cell type classification, complex or simple}
    \item{\code{Area}}{a factor with levels \code{V1} \code{V2} corresponding to the cortical area in which the cell was recorded}
    \item{\code{Cell}}{a factor with levels \code{a}, ..., \code{f} as an identifier of the cell.  The numbers, also, correspond to the panels of Figure 2 from the reference (see below) in which the data were published.}
  }
}
\details{
Extracellular recordings of cells from macaque visual cortical areas V1 and V2.  Simple and complex cells were differentiated by whether or not the response rate, \eqn{f_0}, was greater than the amplitude of modulation at the drift rate, \eqn{f_1}.  If \eqn{f_1 > f_0}, the cell was classified as a simple cell and the value of  \eqn{f_1} is reported; otherwise the cell is classified as a complex cell and the value of \eqn{f_0} is reported.  In both cases, the units are in impulses/second.
}
\source{Data collected in Peter Lennie's lab, analysis 
by Jon Peirce (see reference below)}
\references{
Peirce, J. W. (2007). The potential importance of saturating and supersaturating contrast response functions in visual cortex. \emph{Journal of Vision}, \bold{7(6):13}, 1--10, \url{http://journalofvision.org/7/6/13/}.
}
\examples{
data(CorticalCells)
lattice:::xyplot(Response ~ Contrast | Cell, CorticalCells, 
	type = "b")
}
\keyword{datasets}
