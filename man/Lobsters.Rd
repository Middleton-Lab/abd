% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{Lobsters}
\alias{Lobsters}
\title{Lobster Orientation}
\format{
A data frame with 15 observations of one variable. \describe{
\item{orientation}{} }
}
\source{
Boles, L.C. and K.J. Lohmann. 2003. True navigation and magnetic
maps in spiny lobsters. \emph{Nature} 421: 60-63.
}
\usage{
Lobsters
}
\description{
Orientation of 15 lobsters relative to initial position.
}
\examples{

histogram(~ orientation, Lobsters)
dotplot(~ orientation, Lobsters)

}
\references{
\url{http://www.unc.edu/depts/geomag/PDFGeomag/BolesandLohmann2003.pdf}
}
\keyword{datasets}
