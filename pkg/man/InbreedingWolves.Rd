\name{InbreedingWolves}
\alias{InbreedingWolves}
\docType{data}

\title{
Inbreeding in Wolves
}

\description{
Inbreeding coefficient and the number of pups produced in 24 mated pairs of wolves (\emph{Canis lupus}) from 1983-2002.
}

\usage{data(InbreedingWolves)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{\code{inbreeding.coefficient}}{a numeric vector}
    \item{\code{pups}}{a numeric vector}
  }

}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Liberg, O.H., H. Andrén, H.-C. Pedersen, H. Sand, D. Sejberg, P. Wabakken, M. Åkesson, and S. Bensch. 2005. Severe inbreeding depression in a wild wolf (\emph{Canis lupus}) population. \emph{Biology Letters} 1: 17-20.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(InbreedingWolves)
InbreedingWolves

# Plot with jitter() to separate integer numbers of pups on y axis
plot(jitter(pups) ~ inbreeding.coefficient, data = InbreedingWolves,
  xlab = "Inbreeding Coefficient",
  ylab = "Number of Pups",
  pch = 16, col = "red")

(sum.products <- sum_of_products(
  InbreedingWolves$inbreeding.coefficient,
  InbreedingWolves$pups))

SS.inbreeding <- sum_of_squares(
  InbreedingWolves$inbreeding.coefficient)
SS.pups <- sum_of_squares(InbreedingWolves$pups)
(r <- sum.products / (sqrt(SS.inbreeding) * sqrt(SS.pups)))

# Testing the null hypothesis of zero correlation
n <- nrow(InbreedingWolves)
(SE.r <- sqrt((1 - r^2) / (n - 2)))
(t.stat <- r / SE.r)
2 * pt(t.stat, df = (n - 2))

# Or using rounded values from p. 440
2 * pt(-3.60, 22)

# With cor.test()
cor.test(InbreedingWolves$inbreeding.coefficient,
  InbreedingWolves$pups)
}

\keyword{datasets}
