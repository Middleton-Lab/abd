\name{Toads}
\alias{Toads}
\docType{data}

\title{
Right-handed Toads
}

\description{
Hypothetical probability of a toad being right-handed
}

\usage{data(Toads)}
\format{
  A data frame with 19 observations on the following 2 variables.
  \describe{
    \item{\code{n.toads}}{a numeric vector}
    \item{\code{prob}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(Toads)
Toads
barplot(Toads$prob,
  ylim = c(0, 0.20),
  names.arg = Toads$n,
  cex.names = 0.75)
}

# Using dbinom()
barplot(dbinom(0:18, 18, prob = 0.5),
  ylim = c(0, 0.20),
  names.arg = 1:18,
  cex.names = 0.75)

# Exact two-tailed P-value for n >= 14 right-handed toads
2 * sum(dbinom(14:18, 18, 0.5))

\keyword{datasets}
