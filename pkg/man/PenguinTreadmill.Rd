\name{PenguinTreadmill}
\alias{PenguinTreadmill}
\docType{data}

\title{
Penguin Heart Rate
}

\description{
Slope of regressions of mass-specific metabolic rate on heart rate for three groups of Macaroni Penguins.
}

\usage{data(PenguinTreadmill)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{\code{group}}{a factor with levels \code{BF}, \code{BM}, and \code{MF}}
    \item{\code{slope}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Green, J. A., P. J. Butler, A. J. Woakes, I. L. Boyd and R. L. Holder. 2001. Heart rate and rate of oxygen consumption of exercising macaroni penguins. \emph{Journal of Experimental Biology} 204: 673–684.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(PenguinTreadmill)
str(PenguinTreadmill)
dotplot(slope~group, PenguinTreadmill)
}

\keyword{datasets}
