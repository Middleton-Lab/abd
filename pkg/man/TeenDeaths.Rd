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
    \item{\code{Cause}}{a character vector}
    \item{\code{No.deaths}}{a numeric vector}
  }
}

\source{
Anderson, R.N. 2001. Deaths: Leading causes for 1999. National vital statistics reports 49(11):1-88. National Center for Health Statistics; Hyattsville, Maryland.}

\examples{
data(TeenDeaths)

op <- par(no.readonly = TRUE)
par(mai = c(2, 0.82, 0.25, 0.42))
str(TeenDeaths)
TeenDeaths
barplot(TeenDeaths$No.deaths, names.arg = TeenDeaths$Cause, las = 3, 
	cex.axis = 0.75, cex.names = 0.75, ylim = c(0, 7000),
	ylab = 'Number of Cases (frequency)',
	col = 'red')
par(op)

}

\keyword{datasets}
