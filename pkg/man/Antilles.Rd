\name{Antilles}
\alias{Antilles}
\docType{data}

\title{
Antilles Bird Immigration Dates
}

\description{
Approximate dates of immigration for 37 species of birds in the Lesser Antilles.
}

\usage{data(Antilles)}
\format{A data frame with 37 observations of one variable.
\describe{
\item{\code{immigration.date}}{approximate immigration date (in millions of years)}
}
}

% \details{
% }

\source{
\emph{inferred from} Ricklefs, R.E. and E. Bermingham. 2001. Nonequilibrium diversity dynamics of the Lesser Antillean avifauna. \emph{Science} 294: 1522-1524.
}

\references{
\url{http://www.sciencemag.org/cgi/content/abstract/sci;294/5546/1522}
}

\examples{
data(Antilles)
histogram(~immigration.date, Antilles,n=15)
densityplot(~immigration.date, Antilles)
}

\keyword{datasets}
