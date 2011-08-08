\name{FingerRatio}
\alias{FingerRatio}
\docType{data}

\title{
2D:4D Finger Ratio
}

\description{
The ratio of the lengths of the index finger to the ring finger in 46 males and the number of CAG repeats for each.
}

\usage{data(FingerRatio)}
\format{
  A data frame with 46 observations on the following 2 variables.
  \describe{
    \item{\code{CAGrepeats}}{a numeric vector of the number of CAG repeats}
    \item{\code{finger.ratio}}{a numeric vector of the ratio of digit 2 to digit 4}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Manning, J.T., P.E. Bundred, D.J. Newton, and B.F. Flanagan. 2003. The second to fourth digit ratio and variation in the androgen receptor gene. \emph{Evolution and Human Behavior} 24: 399-405.
}

\references{
\url{http://en.wikipedia.org/wiki/Digit_ratio}
}

\examples{
data(FingerRatio)
str(FingerRatio)
xyplot(finger.ratio ~ CAGrepeats, FingerRatio,
  xlab = "Number of CAG Repeats",
  ylab = "2D:4D Ratio"
)
}

\keyword{datasets}
