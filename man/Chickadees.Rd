% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{Chickadees}
\alias{Chickadees}
\title{Alarm Calls in Chickadees}
\format{
A data frame with 13 observations on the following 3 variables.
\describe{ \item{species}{a character vector} \item{mass}{a
numeric vector} \item{dees}{a numeric vector} }
}
\source{
Templeton, C.N., E. Greene, and K. Davis. 2005. Allometry of alarm
calls: Black-capped Chickadees encode information about predator size.
\emph{Science} 308: 1934-1937.
}
\usage{
Chickadees
}
\description{
Number of "dee" notes per call in Black-capped Chickadees (\emph{Poecile
atricapilla}) for 13 predator species with differing body masses.
}
\examples{

str(Chickadees)
Chickadees

xyplot(dees ~ mass, data = Chickadees,
   xlab = "Predator body mass (kg)",
   ylab = "'Dees' per call", type=c('p','r')
)

}
\references{
\url{http://www.sciencemag.org/cgi/content/short/308/5730/1934}
}
\keyword{datasets}
