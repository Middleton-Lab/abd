\name{Wolves}
\alias{Wolves}
\docType{data}

\title{
Inbreeding in Wolves
}

\description{
Inbreeding coefficient and the number of pups produced in 24 mated pairs of wolves (\emph{Canis lupus}) from 1983-2002.
}

\usage{data(Wolves)}
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
data(Wolves)
Wolves
xyplot( inbreeding.coefficient ~ jutter(pups, amount=0.15, Wolves) )
}

\keyword{datasets}
