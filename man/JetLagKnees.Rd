\name{JetLagKnees}
\alias{JetLagKnees}
\docType{data}

\title{
Circadian Rhythm Phase Shift
}

\description{
Shift in circadian rhythm (hours; \code{shift}) in three light treatments
(\code{treatment}).
}

\usage{data(JetLagKnees)}
\format{
  A data frame with 22 observations on the following 2 variables.
  \describe{
	\item{\code{treatment}}{a factor with levels \code{control}, \code{eyes},
	and \code{knee}}
    \item{\code{shift}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Wright, K.P., Jr. and C.A. Czeisler 2002. Absence of circadian phase resetting in response to bright light behind the knees. \emph{Science} 297: 571.
}

\references{
\url{http://www.sciencemag.org/cgi/content/full/297/5581/571}
}

\examples{
data(JetLagKnees)
JetLagKnees
str(JetLagKnees)

# since data set is small, no need to summarize
xyplot(shift ~ treatment, data = JetLagKnees)
bwplot(shift ~ treatment, data = JetLagKnees)

}

\keyword{datasets}
