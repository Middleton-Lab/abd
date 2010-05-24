\name{NeanderthalBrainSize}
\alias{NeanderthalBrainSize}
\docType{data}

\title{
Cranial Capacity in Neanderthals and Modern Humans
}

\description{
Brain size (\code{lnbrain}) and body mass (\code{lnmass}) in Neanderthals and early modern humans (\code{species}).
}

\usage{data(NeanderthalBrainSize)}
\format{
  A data frame with 39 observations on the following 3 variables.
  \describe{
    \item{\code{lnmass}}{a numeric vector}
    \item{\code{lnbrain}}{a numeric vector}
    \item{\code{species}}{a factor with levels \code{neanderthal} \code{recent}}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Ruff, C.B., E. Trinkaus, and T.W. Holliday. 1997. Body mass and encephalization in Pleistocene \emph{Homo}. \emph{Nature} 387: 173-176.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(NeanderthalBrainSize)
NeanderthalBrainSize
}

\keyword{datasets}
