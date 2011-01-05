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


# Mean, variance, standard deviation
n <- length(GlidingSnakes$undulation.rate)
sum(GlidingSnakes$undulation.rate)/n
mean(GlidingSnakes$undulation.rate)

# Calculate variance by hand
with(GlidingSnakes, sum( (undulation.rate - mean(undulation.rate))^2 ) / (n-1))
# Calculate variance using built-in function
var(GlidingSnakes$undulation.rate)
# Standard deviation equals the square root of the variance
sd(GlidingSnakes$undulation.rate)
sd(GlidingSnakes$undulation.rate)^2 - var(GlidingSnakes$undulation.rate)

# Coefficient of variation
CV <- sd(GlidingSnakes$undulation.rate) / mean(GlidingSnakes$undulation.rate) 
signif(CV,3)
cv(GlidingSnakes$undulation.rate)
}
\keyword{datasets}
