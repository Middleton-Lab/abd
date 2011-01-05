\name{Dioecy}
\alias{Dioecy}
\docType{data}

\title{
Dioecy vs. Monomorphism in Plants
}

\description{
Number of \code{dioecious} and \code{monomorphic} taxa among pairs of closely
related plants.
}

\usage{data(Dioecy)}
\format{
  A data frame with 28 observations on the following 3 variables.
  \describe{
    \item{\code{dioecious}}{a numeric vector}
    \item{\code{monomorphic}}{a numeric vector}
    \item{\code{taxon.pair}}{identifier for pair}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Heilbuth, J.C. 2000. Lower species richness in dioecious clades. \emph{The American Naturalist} 156: 221-241.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(Dioecy)
xyplot(dioecious ~ monomorphic, Dioecy, alpha=.65, pch=16)
}

\keyword{datasets}
