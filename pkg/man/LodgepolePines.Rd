\name{LodgepolePines}
\alias{LodgepolePines}
\docType{data}

\title{
Lodgepole Pine Cone Masses
}

\description{ Masses of cones of lodgepole pines (\code{conemass}) from 16
different habitat types (\code{habitat}) in western North America.  }

\usage{data(LodgepolePines)}

\format{
  A data frame with 16 observations on the following 4 variables.
  \describe{
    \item{\code{habitat}}{a factor with levels: \code{island absent}, \code{island present}, and \code{mainland present}}
    \item{\code{conemass}}{mass of cone}
	\item{\code{location}}{\code{island} or \code{mainland}}
	\item{\code{squirrels}}{\code{absent} or \code{present}}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Edelaar, P. and C.W. Benkman. 2006. Replicated population divergence caused by localised coevolution? A test of three hypotheses in the Red Crossbill-lodgepole pine system. \emph{Journal of Evolutionary Biology} 19: 1651-1659.
}

\references{
\url{http://en.wikipedia.org/wiki/Lodgepole_pine}

\url{http://en.wikipedia.org/wiki/Red_crossbill}
}

\examples{
data(LodgepolePines)
LodgepolePines
str(LodgepolePines)
xyplot(conemass ~ habitat, LodgepolePines)
}

\keyword{datasets}
