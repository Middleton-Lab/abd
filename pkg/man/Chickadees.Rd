\name{Chickadees}
\alias{Chickadees}
\docType{data}

\title{
Alarm Calls in Chickadees
}

\description{
Number of "dee" notes per call in Black-capped Chickadees (\emph{Poecile atricapilla}) for 13 predator species with differing body masses.
}

\usage{data(Chickadees)}
\format{
  A data frame with 13 observations on the following 3 variables.
  \describe{
  \item{\code{species}}{a character vector}
    \item{\code{mass}}{a numeric vector}
    \item{\code{dees}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Templeton, C.N., E. Greene, and K. Davis. 2005. Allometry of alarm calls: Black-capped Chickadees encode information about predator size. \emph{Science} 308: 1934-1937.
}

\references{
\url{http://www.sciencemag.org/cgi/content/short/308/5730/1934}
}

\examples{
data(Chickadees)
str(Chickadees)
Chickadees

xyplot(dees ~ mass, data = Chickadees,
   xlab = "Predator body mass (kg)",
   ylab = "'Dees' per call", type=c('p','r')
)
}

\keyword{datasets}
