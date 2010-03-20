\name{DesertBirds}
\alias{DesertBirds}
\docType{data}
\title{
Desert bird census data
}

\description{
Census data for desert birds.
}

\usage{data(DesertBirds)}

\format{
  A data frame with 43 observations on the following 2 variables.
  \describe{
    \item{\code{Species}}{a character vector}
    \item{\code{Count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Sauer, J.R., J.E. Hines, and J. Fallon. 2003. The North American breeding bird survey, results and analysis 1966-2002. Version 2003.1. USGS Patuxent Wildlife Research Center, Laurel, MD.
}

\references{
\url{http://www.mbr-pwrc.usgs.gov/bbs/}
}

\examples{
data(DesertBirds)

str(DesertBirds)
DesertBirds

hist(DesertBirds$Count,
  breaks = 12,
  ylab = "Frequency (Number of Species)",
  xlab = "Abundance",
  main = "",
  col = "red")

\dontrun{
# With ggplot2
require(ggplot2)
p <- ggplot(DesertBirds, aes(Count))
p + geom_histogram(binwidth = 40, fill = "red") +
  scale_x_continuous("Abundance") +
  scale_y_continuous("Frequency (Number of Species)")}


# Similar to Fig. 2.1-1
Count.sort <- sort(DesertBirds$Count)
Count.relfreq <- cumsum(Count.sort)/max(cumsum(Count.sort))
plot(Count.sort, Count.relfreq,
  type = "l",
  col = "red",
  xlim = c(0, 700),
  xlab = "Species abundance",
  ylab = "Cumulative relative frequency")

\dontrun{
p <- ggplot(data.frame(Count.sort, Count.relfreq), 
  aes(Count.sort, Count.relfreq))
p + geom_step(direction = "vh") +
  scale_x_continuous("Species abundance") +
  scale_y_continuous("Cumulative relative frequency")
}
}

\keyword{datasets}
