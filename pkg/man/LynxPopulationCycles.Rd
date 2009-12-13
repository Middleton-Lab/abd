\name{LynxPopulationCycles}
\alias{LynxPopulationCycles}
\docType{data}

\title{
LynxPopulationCycles
}

\description{
Number of lynx pelts reported in Canada per year from 1752 to 1819.
}

\usage{data(LynxPopulationCycles)}
\format{
  A data frame with 68 observations on the following 2 variables.
  \describe{
    \item{\code{date}}{a numeric vector}
    \item{\code{no.pelts}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Elton, C. and M. Nicholson. 1942. The ten-year cycle in numbers of the lynx in Canada. \emph{Journal of Animal Ecology} 11: 215-244.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(LynxPopulationCycles)

plot(LynxPopulationCycles$date, LynxPopulationCycles$no.pelts,
  type = 'l',
  xlab = 'Year',
  ylab = 'Lynx fur returns')
points(LynxPopulationCycles$date, LynxPopulationCycles$no.pelts,
  col = 'red',
  pch = 16)

\dontrun{
# Alternate form converting to Date class.
Year <- as.Date(paste('01jan', LynxPopulationCycles$date, sep = ''), '\%d\%b\%Y')
LynxPopulationCycles <- cbind(LynxPopulationCycles, Year)

require(ggplot2)
p <- ggplot(LynxPopulationCycles, aes(Year, no.pelts))
p + geom_line() + 
  geom_point(color = 'red') +
  scale_y_continuous('Lynx fur returns') +
  opts(panel.grid.minor = theme_blank()) +
  opts(panel.grid.major = theme_line(size = 0.25, colour = 'white'))
}

}

\keyword{datasets}
