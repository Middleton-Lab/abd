\docType{data}
\name{TreeSeedlings}
\alias{TreeSeedlings}
\title{Tree Seedlings and Sunflecks}
\format{A data frame with 21 observations on the following 2 variables.
\describe{ \item{fleck.duration}{a numeric vector}
\item{growth}{a numeric vector} }}
\source{
  Leakey, A.D.B., J.D. Scholes, and M.C. Press. 2005.
  Physiological and ecological significance of sunflecks
  for dipterocarp seedlings. \emph{Journal of Experimental
  Botany} 56: 469-482.
}
\description{
  Fleck duration (min) and relative seedling growth rate
  (mm/mm/week) for 21 seedlings of \emph{Shorea leprosula}.
}
\examples{
data(TreeSeedlings)
str(TreeSeedlings)
splom(TreeSeedlings)
}
\references{
  \url{http://jxb.oxfordjournals.org/cgi/content/short/56/411/469}
}
\keyword{datasets}

