\name{meanCI}
\alias{meanCI}
\alias{meanPval}
\alias{propCI}
\alias{propPval}
\title{
Confidence Intervals and P-values for a Mean
}
\description{
Convenience functions that extract a confidence interval or p-value.
}
\usage{
meanCI(...)
meanPval(...)
propCI(...)
propPval(...)
}
\arguments{
  \item{\dots}{ arguments to \code{\link{t.test}} or \code{\link{prop.test}}
}
}
\value{
A \code{numeric}.  As a side-effect, additional information about the interval or hypothesis test
is printed to the screen.
}
\author{
Randall Pruim
}


\seealso{
\code{\link{interval}}, \code{\link{pval}}
}

\examples{
bwplot(extra ~ group, data = sleep)
meanCI(extra ~ group, data = sleep)
meanPval(extra ~ group, data = sleep)
propCI(60,100)
propPval(60,100)
}
\keyword{ univar }
