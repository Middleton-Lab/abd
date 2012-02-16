\docType{data}
\name{Mosquitoes}
\alias{Mosquitoes}
\title{Body Size in Anopheles Mosquitoes}
\format{A data frame with 20 observations on the following 2 variables.
\describe{ \item{weight}{a numeric vector} \item{sex}{a
factor with levels \code{female} and \code{male}} }}
\source{
  Lounibos, L.P., N. Nishimura, J. Conn, and R.
  Lourenco-de-Oliveira. 1995. Life history correlates of
  adult size in the malaria vector \emph{Anopheles
  darlingi}. \emph{Memórias do Instituto Oswaldo Cruz} 90:
  769-774.
}
\description{
  Weights of female and male mosquitos (\emph{Anopheles
  darlingi})
}
\examples{
data(Mosquitoes)
xyplot(weight ~ sex, Mosquitoes)
}
\references{
  \url{http://www.bioline.org.br/request?oc95154}
}
\keyword{datasets}

