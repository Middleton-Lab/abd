\name{WillsDebates}
\alias{WillsDebates}
\docType{data}
\title{
Presidential "Wills"
}

\description{
Number of times a presidential candidate said "will," "shall," or "going to" in
presidential debates from 1960-2004 (years incomplete).
}

\usage{data(WillsDebates)}
\format{
  A data frame with 8 observations on the following 6 variables.
  \describe{
    \item{\code{year}}{year of presidential debate(s)}
    \item{\code{winner}}{winner of the popular vote (may not be winner of election)}
    \item{\code{loser}}{loser of popular vote (may not be loser of election)}
    \item{\code{winner.wills}}{number of times will/shall used by winner
	during debates}
    \item{\code{loser.wills}}{number of times will/shall used by loser
	during debates}
    \item{\code{diff.wills}}{difference between number of times 
	will/shall used by two candidates}
  }
}
\examples{
data(WillsDebates)
}

\keyword{datasets}
