\name{AntillesImmigrationDates}
\alias{AntillesImmigrationDates}
\docType{data}

\title{
Antilles immigration dates
}

\description{
Approximate dates of immigration (\code{immigration.date}) for 37 species of birds in the Lesser Antilles.
}

\usage{data(AntillesImmigrationDates)}
\format{
  A data frame with 37 observations on the following variable.
  \describe{
    \item{\code{immigration.date}}{a numeric vector}
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
data(AntillesImmigrationDates)
AntillesImmigrationDates$immigration.date
}

\keyword{datasets}
