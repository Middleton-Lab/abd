\name{GreatTitMalaria}
\alias{GreatTitMalaria}
\docType{data}

\title{
Malaria in Populations of Great Tit
}

\description{ Two-by-two contingency table of malaria (\emph{Plasmodium})
infection status in control and egg-removal populations of Great Tit
(\emph{Parus major}).  }

\usage{data(GreatTitMalaria)}
\format{
  A data frame with 4 observations on the following 3 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{Control} and \code{Egg removal}}
    \item{\code{response}}{a factor with levels \code{Malaria} and \code{No Malaria}}
    \item{\code{count}}{a numeric vector}
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
data(GreatTitMalaria)

str(GreatTitMalaria)
GreatTitMalaria

# Table 2.3-1
GTM.raw <- expand.dft(GreatTitMalaria, "count")

table(GTM.raw)

if(require(vcd)) {
	mosaic(~treatment + response, GTM.raw)
}
\dontrun{
# Fig. 2.3-1
require(ggplot2)
bar <- ggplot(GreatTitMalaria, 
  aes(x = Treatment, y = count, fill = Response))
bar + geom_bar(stat = "identity", position = "dodge")

# Fig. 2.3-2
bar + geom_bar(stat = "identity", position = "fill")}

}

\keyword{datasets}
