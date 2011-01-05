\name{SexualSelection}
\alias{SexualSelection}
\docType{data}

\title{
Sexual Conflict
}

\description{ Number of species in each of two taxa in closely related taxon
pairings and the difference between the two groups. One taxon has multiple
matings (\code{polyandrous.species}) and one has only single matings
(\code{monandrous.species}).
}

\usage{data(SexualSelection)}
\format{
  A data frame with 25 observations on the following 4 variables.
  \describe{
    \item{\code{polyandrous.species}}{a numeric vector}
    \item{\code{monandrous.species}}{a numeric vector}
    \item{\code{difference}}{a numeric vector}
    \item{\code{taxon.pair}}{identifier}
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

histogram(~ difference, SexualSelection, n = 20)

hist(SexualSelection$difference, breaks = 20)

# Calculate the number of tests and the number of negative tests
(n <- length(SexualSelection$difference))
(n.neg <- sum(SexualSelection$difference < 0))

2 * pbinom(q = n.neg, size = n, prob = 0.5)

# With a binomial test
binom.test(n.neg, n, p = 0.5)
}

\keyword{datasets}
