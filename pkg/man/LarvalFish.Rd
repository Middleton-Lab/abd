\name{LarvalFish}
\alias{LarvalFish}
\docType{data}

\title{
Exploited Larval Fish
}

\description{
Age (\code{age}) and coefficient of variation (\code{cv}) in larval fish from exploited and unexploited species (\code{exploited}).
}

\usage{data(LarvalFish)}
\format{
  A data frame with 28 observations on the following 3 variables.
  \describe{
    \item{\code{age}}{a numeric vector}
    \item{\code{cv}}{a numeric vector}
    \item{\code{exploited}}{a factor with levels \code{no} and \code{yes}}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Hsieh, C.H., C.S. Reiss, J.R. Hunter, J.R. Beddington, R.M. May, and G. Sugihara. 2006. Fishing elevates variability in the abundance of exploited species. \emph{Nature} 443: 859-862.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(LarvalFish)
str(LarvalFish)
xyplot(cv ~ age | exploited, LarvalFish)
xyplot(cv ~ age, groups=exploited, LarvalFish)
}

\keyword{datasets}
