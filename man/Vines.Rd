\name{Vines}
\alias{Vines}
\docType{data}

\title{
Climbing Vines
}

\description{
Number of \code{climbing} and \code{nonclimbing} species within closely related general of plants.
}

\usage{data(Vines)}
\format{
  A data frame with 48 observations on the following 2 variables.
  \describe{
    \item{\code{climbing}}{a numeric vector}
    \item{\code{nonclimbing}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Gianoli, E. 2004. Evolution of a climbing habit promotes diversification in
flowering plants. \emph{Proceedings of the Royal Society of London, Series B, Biological Sciences} 271: 2011-2015.
}

\references{
\url{http://rspb.royalsocietypublishing.org/content/271/1552/2011.full.pdf}
}

\examples{
data(Vines)
xyplot(nonclimbing ~ climbing, Vines, scales=list(log=TRUE))
}

\keyword{datasets}
