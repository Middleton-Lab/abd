\name{NematodeLifespan}
\alias{NematodeLifespan}
\docType{data}

\title{
Effects of Trimethadione on Lifespan in Nematodes
}

\description{
\code{lifespan} of the nematode \emph{Caenorhabditis elegans} in control and three experimental \code{treatment}s of the anticonvulsant drug trimethadione.
}

\usage{data(NematodeLifespan)}
\format{
  A data frame with 200 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels: \code{adult}, \code{larva}, \code{larva+adult}, and \code{water}}
    \item{\code{lifespan}}{a numeric vector of lifespan}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Evason, K., C. Huang, I. Yamben, D.F. Covey, and K. Kornfeld. 2005. Anticonvulsant medications extend worm life-span. \emph{Science} 307: 258-262.
}

\references{
\url{http://www.sciencemag.org/cgi/content/abstract/307/5707/258}
}

\examples{
data(NematodeLifespan)
str(NematodeLifespan)
}

\keyword{datasets}
