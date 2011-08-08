\name{Beetles}
\alias{Beetles}
\docType{data}

\title{
Beetle Wings and Horns
}

\description{
Relative size of the horns and wings in 19 female \emph{Onthophagus sagittarius} beetles.
}

\usage{data(Beetles)}
\format{
  A data frame with 19 observations on the following 2 variables.
  \describe{
    \item{\code{horn.size}}{a numeric vector}
    \item{\code{wing.mass}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Emlen, D.J. 2001. Costs and the diversification of exaggerated animal structures. \emph{Science} 291: 1534-1536.
}

\references{
\url{http://www.scienceonline.org/cgi/content/abstract/291/5508/1534}
}

\examples{
data(Beetles)
str(Beetles)
xyplot(wing.mass ~ horn.size, Beetles)
}

\keyword{datasets}
