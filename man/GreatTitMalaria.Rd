\name{GreatTitMalaria}
\alias{GreatTitMalaria}
\docType{data}

\title{
Malaria in Populations of Great Tit
}

\description{
Two-by-two contingency table of malaria (\emph{Plasmodium})
infection status in control and egg-removal populations of Great Tit
(\emph{Parus major}).
}

\usage{data(GreatTitMalaria)}

\format{
  A data frame with 65 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{Control} and \code{Egg removal}}
    \item{\code{response}}{a ordered factor with levels \code{No Malaria} and \code{Malaria}}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Oppliger, A., P. Christe, and H. Richner. 1996. Clutch size and malaria resistance. \emph{Nature} 381: 565.
}

\references{
\url{http://www.nature.com/nature/journal/v381/n6583/abs/381565a0.html}
}

\examples{
str(GreatTitMalaria)

table(GreatTitMalaria)

if(require(vcd)) {
  mosaic(~treatment + response, GreatTitMalaria)
}

}

\keyword{datasets}
