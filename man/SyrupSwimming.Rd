\name{SyrupSwimming}

\alias{SyrupSwimming}

\docType{data}

\title{
Syrup Swimming
}

\description{
Relative swimming speed (speed in syrup / speed in water) for 18 swimmers.
}

\usage{data(SyrupSwimming)}

\format{
A data frame with 18 observations of one variable.
\describe{
\item{\code{relative.speed}}{ratio of speed in syrup to speed in water}
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Gettelfinger, B. and E. L. Cussler. 2004. Will Humans Swim Faster or Slower in Syrup? \emph{AIChE Journal} 50: 2646-2647.
}

\references{
\url{http://www3.interscience.wiley.com/journal/109665380/issue}
}

\examples{
data(SyrupSwimming)
SyrupSwimming
histogram(~relative.speed, SyrupSwimming)
dotplot(~relative.speed, SyrupSwimming)
}

\keyword{datasets}
