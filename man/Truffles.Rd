\docType{data}
\name{Truffles}
\alias{Truffles}
\title{Truffle Distribution}
\format{A data frame with 5 observations on the following 2 variables.
\describe{ \item{truffles}{a numeric vector} \item{count}{a
numeric vector} }}
\source{
  Waters, J.R., K.S. McKelvey, D.L. Luoam, and C.J. Zabel.
  1997. Truffle production in old-growth and mature fir
  stands in northeastern California. \emph{Forest Ecology
  and Management} 96: 155-166.
}
\description{
  Number of truffles per plot for 288 plots in an old
  growth forest in northeastern California.
}
\examples{
data(Truffles)
Truffles
xyplot(count~truffles, Truffles, type='h', lwd=4)
barchart(count~truffles, Truffles, origin=0, horizontal=FALSE)
}
\references{
  \url{http://www.fs.fed.us/psw/publications/watersj/waters2.PDF}
}
\keyword{datasets}

