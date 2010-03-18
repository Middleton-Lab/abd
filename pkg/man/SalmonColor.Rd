\name{SalmonColor}
\alias{SalmonColor}
\docType{data}

\title{
Pacific Salmon Color
}

\description{
Skin color sockeye and kokanee morphs of the Pacific salmon (\emph{Oncorhynchus nerka}) raised raised in a low cartoenoid environment.
}

\usage{data(SalmonColor)}
\format{
  A data frame with 35 observations on the following 2 variables.
  \describe{
    \item{\code{species}}{a factor with levels \code{kokanee} and \code{sockeye}}
    \item{\code{skin.color}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Craig, J.K. and C. Foote. 2001. Countergradient variation and secondary sexual color: phenotypic convergence promotes genetic divergence in carotenoid use between sympatric anadromous and nonanadromous morphs of sockeye salmon (\emph{Oncorhynchus nerka}). \emph{Evolution} 55: 380-391.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(SalmonColor)
## maybe str(SalmonColor) ; plot(SalmonColor) ...
}

\keyword{datasets}
