\name{TeenDeaths}
\alias{TeenDeaths}
\docType{data}
\title{
Causes of Teenage Deaths
}

\description{
Data from Table 1 (p. 14) on causes of death for all races, both sexes, ages 15-19.
}

\usage{data(TeenDeaths)}

\format{
  A data frame with 11 observations on the following 2 variables.
  \describe{
    \item{\code{cause}}{a character vector}
    \item{\code{deaths}}{a numeric vector}
  }
}

\source{
Anderson, R.N. 2001. Deaths: Leading causes for 1999. \emph{National vital statistics reports} 49(11):1-88. National Center for Health Statistics; Hyattsville, Maryland.}

\examples{
data(TeenDeaths)

str(TeenDeaths)
TeenDeaths
barchart(deaths~cause, TeenDeaths, 
	horizontal=FALSE,
	ylab="Number of Deaths",
	xlab="Cause of Death",
	scales=list(x=list(rot=45)))
barchart(deaths~ordered(cause, levels=cause), TeenDeaths, 
	horizontal=FALSE,
	ylab="Number of Deaths",
	xlab="Cause of Death",
	scales=list(x=list(rot=45))
	)

}
\keyword{datasets}
