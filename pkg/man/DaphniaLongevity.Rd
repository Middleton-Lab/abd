\name{DaphniaLongevity}
\alias{DaphniaLongevity}
\docType{data}

\title{
Daphnia Longevity
}

\description{
Number of spores and host longevity in the crustacean \emph{Daphnia magna}.
}

\usage{data(DaphniaLongevity)}
\format{
  A data frame with 32 observations on the following 2 variables.
  \describe{
    \item{\code{longevity}}{a numeric vector}
    \item{\code{sqrt.spores}}{a numeric vector}
  }
}

%\details{
%}

\source{
Jensen, K.H., T.J. Little, A. Skorping, and D. Ebert. 2006. Empirical support for optimal virulence in a castrating parasite. \emph{PLoS Biology} 4(7): e197
}

\references{
\url{http://www.plosbiology.org/article/info:doi/10.1371/journal.pbio.0040197}
}

\examples{
data(DaphniaLongevity)
str(DaphniaLongevity)
xyplot(sqrt.spores ~ longevity, DaphniaLongevity)
}

\keyword{datasets}
