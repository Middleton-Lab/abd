\docType{data}
\name{Blackbirds}
\alias{Blackbirds}
\title{Testosterone Levels in Blackbirds}
\format{A data frame with 13 observations on the following 6 variables.
\describe{ \item{before}{a numeric vector} \item{after}{a
numeric vector} \item{log.before}{a numeric vector}
\item{log.after}{a numeric vector} \item{diff.in.logs}{a
numeric vector} \item{diff}{a numeric vector} }}
\source{
  Hasselquist, D., J.A. Marsh, P.W. Sherman, and J.C.
  Wingfield. 1999. Is avian immunocompetence suppressed by
  testosterone? \emph{Behavioral Ecology and Sociobiology}
  45: 167-175.
}
\description{
  Experimental manipulation of testosterone levels in male
  Red-winged Blackbirds (\emph{Agelaius phoeniceus}) and
  resulting changes in antibody levels
}
\examples{
data(Blackbirds)
Blackbirds
xyplot(log.after ~ log.before, data = Blackbirds,
  ylab = "log Antibody production after implant",
  xlab = "log Antibody production before implant"
)
}
\keyword{datasets}

