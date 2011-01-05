\name{Tobacco2}
\alias{Tobacco2}
\docType{data}
\title{
Flower Length in Tobacco Plants
}
\description{
Distribution of flow lengths in F1 and F2 populations of \emph{Nicotiana}.
}

\usage{data(Tobacco2)}
\format{
  A data frame with 617 observations on the following 2 variables.
  \describe{
    \item{\code{flower.length}}{a numeric vector}
    \item{\code{generation}}{a factor with levels \code{F1} \code{F2}}
  }
}
%\details{
%%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
East, E.M. 1916. Studies On Size Inheritance in \emph{Nicotiana}. \emph{Genetics} 1: 164-176.
}

\references{
\url{http://www.genetics.org/content/vol1/issue2/}

\url{http://en.wikipedia.org/wiki/Nicotiana}
}
\seealso{\code{\link{Tobacco}}}

\examples{
data(Tobacco2)
xtabs(~flower.length + generation, Tobacco2)
bwplot(flower.length ~ generation, Tobacco2)
}

\keyword{datasets}
