\docType{data}
\name{Zooplankton}
\alias{Zooplankton}
\title{Zooplankton Depredation}
\format{A data frame with 15 observations on the following 3 variables.
\describe{ \item{treatment}{a factor with levels \code{control},
\code{high}, and \code{low}} \item{zooplankton}{a numeric vector}
\item{block}{a numeric vector} }}
\source{
  \emph{inferred from} Svanbäck, R. and D.I. Bolnick. 2007.
  Intraspecific competition drives increased resource use
  diversity within a natural population. \emph{Proceedings
  of the Royal Society of London Series B, Biological
  Sciences} 274: 839-844.
}
\description{
  Diversity of zooplankton (\code{zooplankton}) prey in
  each of 5 replicate blocks (\code{block}) of three
  treatment levels (\code{treatment}). By default,
  \code{block} is not coded as a factor.
}
\examples{
data(Zooplankton)
Zooplankton

Zooplankton$block <- factor(Zooplankton$block)
str(Zooplankton)

aov.fit <- aov(zooplankton ~ block + treatment,
  data = Zooplankton)
summary(aov.fit)
}
\keyword{datasets}

