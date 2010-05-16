\name{SexualSelection}
\alias{SexualSelection}
\docType{data}

\title{
Sexual Conflict
}

\description{
Number of species in each of two taxa in closely related taxon pairings and the difference between the two groups. One taxon has multiple matings (\code{Polyandrous.clade.richness}) and one has only single matings (\code{Monandrous.clade.richness}).
}

\usage{data(SexualSelection)}
\format{
  A data frame with 25 observations on the following 3 variables.
  \describe{
    \item{\code{Polyandrous.clade.richness}}{a numeric vector}
    \item{\code{Monandrous.clade.richness}}{a numeric vector}
    \item{\code{Difference}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Arnqvist, G., M. Edvardsson, U. Friberg, and T. Nilsson. 2000. Sexual conflict promotes speciation in insects. \emph{Proceedings of the National Academy of Sciences (USA)} 97: 10460-10464.
}

\references{
\url{http://www.egs.uu.se/evbiol/Persons/Urban/1.pdf}
}

\examples{
data(SexualSelection)
SexualSelection

hist(SexualSelection$Difference, breaks = 20)

# Calculate the number of tests and the number of negative tests
(n <- length(SexualSelection$Difference))
(n.neg <- sum(SexualSelection$Difference < 0))

2 * pbinom(q = n.neg, size = n, prob = 0.5)

# With a binomial test
binom.test(n.neg, n, p = 0.5)
}

\keyword{datasets}
