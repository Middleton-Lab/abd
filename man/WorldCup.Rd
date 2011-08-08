\name{WorldCup}
\alias{WorldCup}
\docType{data}

\title{
World Cup Goals
}

\description{
Number of goals per team during the 2002 World Cup.
}

\usage{data(WorldCup)}
\format{
  A data frame with 7 observations on the following 2 variables.
  \describe{
    \item{\code{score}}{a numeric vector}
    \item{\code{count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(WorldCup)
xyplot(count ~ score, WorldCup, type='h', lwd=4)
barchart(count ~ score, WorldCup, horizontal=FALSE)
}

\keyword{datasets}
