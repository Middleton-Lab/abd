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
propCI(..., method=c('binom.test','prop.test', 'wilson', 'agresti-coull'))
propPval(..., method=c('binom.test','prop.test'))
}
\arguments{
  \item{\dots}{ arguments to \code{\link{t.test}} or \code{\link{prop.test}}}
  \item{method}{ an \code{integer} or \code{character} naming which method to use }
}
\value{
A \code{numeric}.  As a side-effect, additional information about the interval
or hypothesis test is printed to the screen.
}
\author{
Randall Pruim (\email{rpruim@calvin.edu})
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
