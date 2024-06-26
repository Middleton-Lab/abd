% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{SpiderSpeed}
\alias{SpiderSpeed}
\title{Spider Running Speeds after Amputation}
\format{
A data frame with 32 observations on the following 2 variables.
\describe{ \item{speed.before}{speed (cm/s) before amputation }
\item{speed.after}{speed (cm/s) after amputation } }
}
\source{
Ramos, M., D.J. Irschick, and T.E. Christenson. 2004. Overcoming an
evolutionary conflict: Removal of a reproductive organ greatly increases
locomotor performance. \emph{Proceedings of the National Academy of Sciences
(USA)} 101: 4883-4887.
}
\usage{
SpiderSpeed
}
\description{
Data on speed before and after amputation of a pedipalp in the spider genus
\emph{Tidarren}.
}
\examples{
xyplot(speed.after ~ speed.before, SpiderSpeed)
favstats(SpiderSpeed$speed.before)
favstats(SpiderSpeed$speed.after)
favstats(SpiderSpeed$speed.after - SpiderSpeed$speed.before)

}
\references{
\url{http://en.wikipedia.org/wiki/Pedipalp}, 
\url{http://en.wikipedia.org/wiki/Tidarren}, 
\url{http://www.pnas.org/content/101/14/4883}
}
\keyword{datasets}
