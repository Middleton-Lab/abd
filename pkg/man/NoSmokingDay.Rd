\name{NoSmokingDay}
\alias{NoSmokingDay}
\docType{data}

\title{
No Smoking Day
}

\description{
Number of workplace injuries on No Smoking Day (\code{Injuries.on.NSD}) compared to the same Wednesday in the previous year (\code{Injuries.before.NSD}) for 1987-1996.
}

\usage{data(NoSmokingDay)}
\format{
  A data frame with 10 observations on the following 3 variables.
  \describe{
    \item{\code{year}}{a numeric vector}
    \item{\code{Injuries.before.NSD}}{a numeric vector}
    \item{\code{Injuries.on.NSD}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Waters, A.J., M.J. Jarvis, and S.R. Sutton. 1998. Nicotine withdrawal and accident rates. \emph{Nature} 394: 137.
}

\references{
\url{http://www.nosmokingday.org.uk/}
}

\examples{
data(NoSmokingDay)
NoSmokingDay
}

\keyword{datasets}
