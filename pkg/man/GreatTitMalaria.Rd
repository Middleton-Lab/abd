\name{GreatTitMalaria}
\alias{GreatTitMalaria}
\docType{data}

\title{
GreatTitMalaria
}

\description{
Two-by-two contingency table of malaria (\emph{Plasmodium}) infection frequencies in control and egg-removal populations of Great Tit (\emph{Parus major}).
}

\usage{data(GreatTitMalaria)}
\format{
  A data frame with 4 observations on the following 3 variables.
  \describe{
    \item{\code{Treatment}}{a factor with levels \code{Control} and \code{Egg.removal}}
    \item{\code{Response}}{a factor with levels \code{Malaria} and \code{No.Malaria}}
    \item{\code{Frequency}}{a numeric vector}
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
# see https://stat.ethz.ch/pipermail/r-help/2009-January/185561.html
# for discussion of expand.dft(). Modified for GreatTitMalaria data.
expand.dft <- function(x, na.strings = 'NA', as.is = FALSE, dec = '.'){
  DF <- sapply(1:nrow(x), function(i) x[rep(i, each = x$Freq[i]), ], 
    simplify = FALSE)
  DF <- subset(do.call('rbind', DF), select = -Frequency)
  for (i in 1:ncol(DF)){
    DF[[i]] <- type.convert(as.character(DF[[i]]),
      na.strings = na.strings, as.is = as.is, dec = dec)
  }
DF
} 

GTM.raw <- expand.dft(GreatTitMalaria)

require(gmodels)
CrossTable(GTM.raw$Treatment, GTM.raw$Response,
  expected = FALSE,
  prop.r = FALSE,
  prop.c = FALSE,
  prop.chisq = FALSE, 
  prop.t = FALSE)

# Fig. 2.3-1
\dontrun{
require(ggplot2)
bar <- ggplot(GreatTitMalaria, aes(x = Treatment, y = Frequency, fill = Response))
bar + geom_bar(stat = 'identity', position = 'dodge')

# Fig. 2.3-2
bar + geom_bar(stat = 'identity', position = 'fill')
}

}

\keyword{datasets}
