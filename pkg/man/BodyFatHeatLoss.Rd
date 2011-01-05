\name{BodyFatHeatLoss}
\alias{BodyFatHeatLoss}
\docType{data}

\title{
Heat Loss and Body Fat
}

\description{
Heat loss during exercise and relative body fat in 12 boys.
}

\usage{data(BodyFatHeatLoss)}
\format{
  A data frame with 12 observations on the following 2 variables.
  \describe{
    \item{\code{leanness}}{a numeric vector}
    \item{\code{lossrate}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Sloan, R.E.G. and W.R. Keatinge. 1973. Cooling rates of young people swimming in cold water. \emph{Journal of Applied Physiology} 35: 371-375.
}

\references{
\url{http://www.ncbi.nlm.nih.gov/pubmed/4732330}
}

\examples{
data(BodyFatHeatLoss)
xyplot(lossrate ~ leanness, BodyFatHeatLoss)
}

\keyword{datasets}
