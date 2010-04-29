\name{GreenSpaceBiodiversity}
\alias{GreenSpaceBiodiversity}
\docType{data}

\title{
Diversity in Urban Green Space
}

\description{
Measures of biodiversity in 15 urban green spaces in Sheffield, England.
}

\usage{data(GreenSpaceBiodiversity)}
\format{
  A data frame with 15 observations on the following 6 variables.
  \describe{
    \item{\code{site}}{a factor with levels \code{A} - \code{O}}
    \item{\code{attachment}}{a numeric vector}
    \item{\code{area.ha}}{a numeric vector}
    \item{\code{butterfly}}{a numeric vector}
    \item{\code{bird}}{a numeric vector}
    \item{\code{ln.plant}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Fuller, R.A., K.N. Irvine, P. Devine-Wright, P.H. Warren, and K.J. Gaston. 2007. Psychological benefits of greenspace increase with biodiversity. \emph{Biology Letters} 3: 390-394.
}

\references{
\url{http://rsbl.royalsocietypublishing.org/content/3/4/390.abstract}
}

\examples{
data(GreenSpaceBiodiversity)
str(GreenSpaceBiodiversity)
GreenSpaceBiodiversity
}

\keyword{datasets}
