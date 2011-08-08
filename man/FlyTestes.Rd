\name{FlyTestes}
\alias{FlyTestes}
\docType{data}

\title{
Testes Size in Flies
}

\description{
Testes size (square mm; \code{Testes.area}) in 8 populations of common yellow dung flies (\emph{Scathophaga stercoraria}) with different mating systems (\code{Mating.system}).
}

\usage{data(FlyTestes)}

\format{
  A data frame with 8 observations on the following 2 variables.
  \describe{
    \item{\code{mating}}{a factor with levels \code{Monogamous} \code{Polyandrous}}
    \item{\code{testes.area}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Hosken, D.J. and P.I. Ward. 2001.  Experimental evidence for testis size evolution via sperm competition. \emph{Ecology Letters} 4: 10-13.
}

\references{
\url{http://en.wikipedia.org/wiki/Scathophaga_stercoraria}
}

\examples{
data(FlyTestes)
str(FlyTestes)
FlyTestes
}

\keyword{datasets}
