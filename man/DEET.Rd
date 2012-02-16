\docType{data}
\name{DEET}
\alias{DEET}
\title{DEET and Mosquito Bites}
\format{A data frame with 52 observations on the following 2 variables.
\describe{ \item{dose}{a numeric vector} \item{bites}{a
numeric vector} }}
\source{
  Golenda, C.F., V.B. Solberg, R. Burge, J.M. Gambel, and
  R.A. Wirtz. 1999. Gender-related efficacy difference to
  an extended duration formulation of topical
  N,N-diethyl-\emph{m}-toluamide (DEET). \emph{American
  Journal of Tropical Medicine and Hygiene} 60: 654-657.
}
\description{
  Administered dose of DEET and number of mosquito bites
  for 52 women.
}
\examples{
data(DEET)
str(DEET)
xyplot(bites ~ dose, DEET)
}
\keyword{datasets}

