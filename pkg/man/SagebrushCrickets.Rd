\name{SagebrushCrickets}
\alias{SagebrushCrickets}
\docType{data}

\title{
Sagebrush Cricket Mating Times
}

\description{
Time to mating (\code{time.to.mating}) in fed and unfed (\code{treatment}) 
sagebrush crickets (\emph{Cyphoderris strepitans}).
}

\usage{data(SagebrushCrickets)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels: \code{fed} and \code{starved}}
    \item{\code{time.to.mating}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Chadwick Johnson, J., T.M. Ivy, and S.K. Sakaluk. 1999. Female remating propensity contingent on sexual cannibalism in sagebrush crickets, \emph{Cyphoderris strepitans}: a mechanism of cryptic female choice. \emph{Behavioral Ecology} 10: 227-233.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(SagebrushCrickets)
SagebrushCrickets
str(SagebrushCrickets)
}

\keyword{datasets}
