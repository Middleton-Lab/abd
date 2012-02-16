\docType{data}
\name{Lefthanded}
\alias{Lefthanded}
\title{Left-handedness and Rates of Violence}
\format{A data frame with 8 observations on the following 2 variables.
\describe{ \item{percent.left}{a numeric vector}
\item{murder.rate}{a numeric vector} }}
\source{
  Faurie, C. and M. Raymond. 2005. Handedness, homicide and
  negative frequency-dependent selection. \emph{Proceedings
  of the Royal Society of London B} 272: 25-28.
}
\description{
  Prevalence of left-handedness (\code{percent.left}) and
  homicide rates (\code{murder}) for 8 societies.
}
\examples{
data(Lefthanded)
str(Lefthanded)
xyplot(murder.rate ~ percent.left, Lefthanded)
}
\references{
  \url{http://rspb.royalsocietypublishing.org/content/272/1558/25.abstract}
}
\keyword{datasets}

