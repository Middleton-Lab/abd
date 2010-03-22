\name{Newts}
\alias{Newts}
\docType{data}

\title{
Tetrodotoxin Resistance in Garter Snakes
}

\description{
Percent reduction in crawl speed (\code{whole.animal.resistance}) in the garter snake after injection of the neurotoxin tetrodotoxin from the rough-skinned newt (\emph{Taricha granulosa}).
}

\usage{data(Newts)}
\format{
  A data frame with 12 observations on the following 2 variables.
  \describe{
    \item{\code{locality}}{a factor with levels: \code{Benton} and \code{Warrenton}}
    \item{\code{whole.animal.resistance}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Geffeney, S., E.D. Brodie, Jr., P.C. Ruben, and E.D. Brodie III. 2002. Mechanisms of adaptation in a predator-prey arms race: TTX-resistant sodium channels. \emph{Science} 297: 1336-1339.
}

\references{
\url{http://www.sciencemag.org/cgi/content/abstract/297/5585/1336}
}

\examples{
data(Newts)
Newts
}

\keyword{datasets}
