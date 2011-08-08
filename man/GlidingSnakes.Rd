\name{GlidingSnakes}
\alias{GlidingSnakes}
\docType{data}

\title{
GlidingSnakes
}

\description{ Undulation rate (\emph{Hz}) of 8 paradise tree snakes
(\emph{Chrysopelea paradisi}).  }

\usage{data(GlidingSnakes)}
\format{
A data frame with eight observations of one variable.
\describe{
\item{\code{undulation.rate}}{undulation rate}
}
}

\source{
Socha, J.J. 2002. Gliding flight in the paradise tree snake. \emph{Nature} 418: 603-604.
}

\references{
\url{http://www.nature.com/nature/journal/v418/n6898/abs/418603a.html}

\url{http://www.flyingsnake.org/}
}

\examples{
data(GlidingSnakes)

histogram(~undulation.rate , data=GlidingSnakes, n=7,
  xlab = "Undulation rate (Hz)",
  type='count')

}
\keyword{datasets}
