\name{favstats}
\alias{favstats}
\title{Some favorite statistical summaries }
\description{
Computes median, IQR, mean, standard deviation, and variance of a data vector.
}
\usage{
favstats(x, na.rm = TRUE)
}
\arguments{
  \item{x}{ numeric vector }
  \item{na.rm}{ boolean indicating whether missing data should be ignored }
}
\value{A vector of (named) statistical summaries}

\examples{
favstats(1:10)
favstats(faithful$eruptions)
}
\keyword{ stats }
\keyword{ univar }
