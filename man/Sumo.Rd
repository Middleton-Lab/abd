\name{Sumo}
\alias{Sumo}
\docType{data}

\title{
Sumo Wrestling Wins
}

\description{
Counts of number of wins for sumo wrestlers.
}

\usage{data(Sumo)}
\format{
  A data frame with 16 observations on the following 2 variables.
  \describe{
    \item{\code{wins}}{a numeric vector}
    \item{\code{count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Duggan, M. and S.D. Leavitt. 2002. Winning isn't everything: Corruption in sumo wrestling. \emph{The American Economic Review} 92: 1594-1605.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(Sumo)
xyplot(count ~ wins, Sumo, type='h', lwd=4)
}

\keyword{datasets}
