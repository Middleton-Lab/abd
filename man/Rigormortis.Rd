\docType{data}
\name{Rigormortis}
\alias{Rigormortis}
\title{Rigormortis and Time of Death}
\format{A data frame with 12 observations on the following 2 variables.
\describe{ \item{hours}{a numeric vector} \item{count}{a
numeric vector} }}
\source{
  Pounder, D.J. 1995. Postmortem changes and time of death.
  University of Dundee.
  \url{http://www.dundee.ac.uk/forensicmedicine/llb/timedeath.htm}.
}
\description{
  Number of bodies reaching rigormortis in each hour after
  death.
}
\examples{
data(Rigormortis)
xyplot(count~hours, Rigormortis, type='h', lwd=3)
barchart(count ~ hours, Rigormortis, horizontal=FALSE, origin=0)
}
\keyword{datasets}

