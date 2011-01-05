\name{EarwigForceps}
\alias{EarwigForceps}
\docType{data}

\title{
Earwig Density and Forceps
}

\description{
Earwig (\emph{Forficula auricularia}) density and the proportion of trapped earwigs with abdominal forceps (used for fighting and courtship).
}

\usage{data(EarwigForceps)}
\format{
  A data frame with 7 observations on the following 2 variables.
  \describe{
    \item{\code{density}}{a numeric vector}
    \item{\code{proportion.forceps}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Tomkins, J.L.  and G.S. Brown. 2004. Population density drives the local evolution of a threshold dimorphism. \emph{Nature} 431: 1099-1103.
}

\references{
\url{http://en.wikipedia.org/wiki/Forficula_auricularia}

\url{http://www.arkive.org/common-european-earwig/forficula-auricularia/}

\url{http://eol.org/pages/473785}
}

\examples{
data(EarwigForceps)
xyplot(proportion.forceps ~ density, data=EarwigForceps, type='h', lwd=6)

}

\keyword{datasets}
