\docType{data}
\name{VasopressinVoles}
\alias{VasopressinVoles}
\title{Vasopressin Manipulation in the Meadow Vole}
\format{A data frame with 31 observations on the following 2 variables.
\describe{ \item{treatment}{a factor with levels \code{control} and
\code{enhanced}} \item{percent}{a numeric vector} }}
\source{
  \emph{inferred from} Lim, M.M., Z. Wang, D.E. Olazabal,
  X. Ren, E.F. Terwilliger, and L.J. Young. 2004. Enhanced
  partner preference in a promiscuous species by
  manipulating the expression of a single gene.
  \emph{Nature} 429: 754-757.
}
\description{
  Time spent with a female (\code{percent}) in control and
  vasopressin-enhanced groups (\code{treatment}) of meadow
  voles (\emph{Microtus pennsylvanicus}).
}
\examples{
data(VasopressinVoles)
xyplot(percent ~ treatment, VasopressinVoles, type=c('p','a'))
bwplot(percent ~ treatment, VasopressinVoles)
}
\references{
  \url{http://www.nature.com/nature/journal/v429/n6993/full/nature02539.html}
}
\keyword{datasets}

