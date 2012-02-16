\docType{data}
\name{AlgaeCO2}
\alias{AlgaeCO2}
\title{Carbon Dioxide and Growth Rate in Algae}
\format{A data frame with 14 observations on the following 2 variables.
\describe{ \item{treatment}{a factor with levels \code{Normal} and
\code{High}} \item{growthrate}{a numeric vector} }}
\source{
  Collins, S. and G. Bell. 2004. Phenotypic consequences of
  1,000 generations of selection at elevated
  CO\eqn{_{2}}{2} in a green alga. \emph{Nature} 431:
  566-569.
}
\description{
  Growth rates of the unicellular alga \emph{Chlamydomonas}
  after 1,000 generations of selection under \code{High}
  and \code{Normal} levels of carbon dioxide.
}
\examples{
data(AlgaeCO2)
AlgaeCO2
xyplot(growthrate ~ treatment, AlgaeCO2, type = c('p', 'a'))
}
\references{
  \url{http://www.nature.com/nature/journal/v431/n7008/full/nature02945.html}
}
\keyword{datasets}

