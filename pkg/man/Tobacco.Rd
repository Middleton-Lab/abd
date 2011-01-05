\name{Tobacco}
\alias{Tobacco}
\docType{data}

\title{
Flower Length in Tobacco Plants
}

\description{
Distribution of flow lengths in F1 and F2 populations of \emph{Nicotiana}.
}

\usage{data(Tobacco)}
\format{
  A data frame with 13 observations on the following 3 variables.
  \describe{
    \item{\code{flower.length}}{a numeric vector of flower length in mm}
    \item{\code{f1.count}}{a numeric vector of the number of F1 plants with flower lengths in this size range}
    \item{\code{f2.count}}{a numeric vector of the number of F2 plants with flower lengths in this size range}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
East, E.M. 1916. Studies On Size Inheritance in \emph{Nicotiana}. \emph{Genetics} 1: 164-176.
}

\references{
\url{http://www.genetics.org/content/vol1/issue2/}

\url{http://en.wikipedia.org/wiki/Nicotiana}
}

\seealso{\code{\link{Tobacco2}}}
\examples{
data(Tobacco)
Tobacco
}

\keyword{datasets}
