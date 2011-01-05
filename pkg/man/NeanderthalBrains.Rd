\name{NeanderthalBrains}
\alias{NeanderthalBrains}
\docType{data}

\title{
Cranial Capacity in Neanderthals and Modern Humans
}

\description{
Brain size (\code{lnbrain}) and body mass (\code{lnmass}) in Neanderthals and early modern humans (\code{species}).
}

\usage{data(NeanderthalBrains)}
\format{
  A data frame with 39 observations on the following 3 variables.
  \describe{
    \item{\code{ln.mass}}{log of body mass (kg)}
    \item{\code{ln.brain}}{log of brain size}
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
data(NeanderthalBrains)
xyplot(ln.brain ~ ln.mass, data=NeanderthalBrains, groups=species)
}

\keyword{datasets}
