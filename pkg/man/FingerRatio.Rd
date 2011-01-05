\name{FingerRatio}
\alias{FingerRatio}
\docType{data}

\title{
2D:4D Finger Ratio
}

\description{
The ratio of the lengths of the index finger to the ring finger in 46 males and the number of CAG repeats for each.
}

\usage{data(FingerRatio)}
\format{
  A data frame with 46 observations on the following 2 variables.
  \describe{
    \item{\code{CAGrepeats}}{a numeric vector of the number of CAG repeats}
    \item{\code{finger.ratio}}{a numeric vector of the ratio of digit 2 to digit 4}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Manning, J.T., P.E. Bundred, D.J. Newton, and B.F. Flanagan. 2003. The second to fourth digit ratio and variation in the androgen receptor gene. \emph{Evolution and Human Behavior} 24: 399-405.
}

\references{
\url{http://en.wikipedia.org/wiki/Digit_ratio}
}

\examples{
data(FingerRatio)
str(FingerRatio)

plot(FingerRatio$CAGrepeats,
  FingerRatio$finger.ratio,
  xlab = "Number of CAG Repeats",
  ylab = "2D:4D Ratio",
  pch = 16, col = "red")

# Shorten the names a bit
repeats <- FingerRatio$CAGrepeats
ratio <- FingerRatio$finger.ratio

(sum.products <- sum((repeats - mean(repeats)) *
  (ratio - mean(ratio))))
(SS.repeats <- sum((repeats - mean(repeats))^2))
(SS.ratio <- sum((ratio - mean(ratio))^2))

(r <- sum.products / (sqrt(SS.repeats) * sqrt(SS.ratio)))

# Functions from abd package to calculate the same values
sum.products <- sum_of_products(repeats, ratio)
SS.repeats <- sum_of_squares(repeats)
SS.ratio <- sum_of_squares(ratio)
sum.products / (sqrt(SS.repeats) * sqrt(SS.ratio))

# cor() does the calculation in one step.
# Default is Pearson's correlation.
cor(FingerRatio$CAGrepeats,
  FingerRatio$finger.ratio)

# Standard error of r.
# Use nrow() to get the number of observations.
n <- nrow(FingerRatio)
(SE.r <- sqrt((1 - r^2) / (n - 2)))

# Approximate confidence interval
z <- 0.5 * log((1 + r) / (1 - r))
s.z <- sqrt(1 / (n - 3))
z.crit <- qnorm((1 - 0.05/2))

# Lower and upper 95% CIs
(ci.lower <- z - z.crit * s.z)
(ci.upper <- z + z.crit * s.z)

# Backtransformation
(exp(2 * ci.lower) - 1) / (exp(2 * ci.lower) + 1)
(exp(2 * ci.upper) - 1) / (exp(2 * ci.upper) + 1)
}

\keyword{datasets}
