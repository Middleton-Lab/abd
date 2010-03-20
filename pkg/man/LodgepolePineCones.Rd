\name{LodgepolePineCones}
\alias{LodgepolePineCones}
\docType{data}

\title{
Lodgepole Pine Cone Masses
}

\description{
Masses of cones of lodgepole pines (\code{conemass}) from 16 different habitat types (\code{habitat}) in western North America.
}

\usage{data(LodgepolePineCones)}

\format{
  A data frame with 16 observations on the following 2 variables.
  \describe{
    \item{\code{habitat}}{a factor with levels: \code{island.absent}, \code{island.present}, and \code{mainland.present}}
    \item{\code{conemass}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Edelaar, P. and C. W. Benkman. 2006. Replicated population divergence caused by localised coevolution? A test of three hypotheses in the red crossbill-lodgepole pine system. \emph{Journal of Evolutionary Biology} 19: 1651-1659.
}

%\references{
%\url{http://www.uwyo.edu/benkman/PDFs%2520of%2520papers/Edelaar&Benkman2006.pdf}
%}

\examples{
data(LodgepolePineCones)
LodgepolePineCones
str(LodgepolePineCones)
}

\keyword{datasets}
