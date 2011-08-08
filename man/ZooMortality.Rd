\name{ZooMortality}
\alias{ZooMortality}
\docType{data}

\title{
Home Range Size and Mortality
}
   
\description{
Home range size (\eqn{\log_{10}}{log10}) and captive infant mortality (\%) for 20 species of carnivores.
}

\usage{data(ZooMortality)}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{log.homerange}}{a numeric vector}
    \item{\code{mortality}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Clubb, R. and G. Mason. 2003. Captivity effects on wide ranging carnivores. \emph{Nature} 425: 473-474.
}

\references{
\url{http://www.nature.com/nature/journal/v425/n6957/full/425473a.html}
}

\examples{
data(ZooMortality)
str(ZooMortality)
}

\keyword{datasets}
