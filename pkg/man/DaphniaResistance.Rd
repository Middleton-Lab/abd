\name{DaphniaResistance}
\alias{DaphniaResistance}
\docType{data}

\title{
Daphnia Resistance to Cyanobacteria
}

\description{
Resistance of \emph{Daphnia} eggs to different levels of cyanobacteria (\code{cyandensity}) from 1962-1997.
}

\usage{data(DaphniaResistance)}
\format{
  A data frame with 32 observations on the following 2 variables.
  \describe{
    \item{\code{density}}{a factor with levels: \code{high}, \code{low}, and \code{med}}
    \item{\code{resistance}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Hairston, N.G., Jr., W. Lampert, C.E. Cáceres, C.L. Holtmeier, L.J. Weider, U. Gaedke, J.M. Fischer, J.A. Fox, and D.M. Post. 1999. Dormant eggs record rapid evolution. \emph{Nature} 401: 446.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(DaphniaResistance)
str(DaphniaResistance)


bwplot(resistance ~ density, DaphniaResistance)
# with such a small data set, we can display all the data rather than a summary
xyplot(resistance ~ density, DaphniaResistance)
histogram( ~ resistance| density, DaphniaResistance, 
	strip=FALSE, strip.left = TRUE,
	layout=c(1,3)
	)
}

\keyword{datasets}
