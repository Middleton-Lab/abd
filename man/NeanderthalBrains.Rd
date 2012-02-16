\docType{data}
\name{NeanderthalBrains}
\alias{NeanderthalBrains}
\title{Cranial Capacity in Neanderthals and Modern Humans}
\format{A data frame with 39 observations on the following 3 variables.
\describe{ \item{ln.mass}{log of body mass (kg)}
\item{ln.brain}{log of brain size} \item{species}{a factor
with levels \code{neanderthal} \code{recent}} }}
\source{
  Ruff, C.B., E. Trinkaus, and T.W. Holliday. 1997. Body
  mass and encephalization in Pleistocene \emph{Homo}.
  \emph{Nature} 387: 173-176.
}
\description{
  Brain size (\code{lnbrain}) and body mass (\code{lnmass})
  in Neanderthals and early modern humans (\code{species}).
}
\examples{
data(NeanderthalBrains)
xyplot(ln.brain ~ ln.mass, data=NeanderthalBrains, groups=species)
}
\keyword{datasets}

