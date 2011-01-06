\name{Lynx}
\alias{Lynx}
\docType{data}

\title{
Population Cycles of Lynx in Canada 1752-1819
}

\description{
Number of lynx pelts (\code{pelts}) reported in Canada per year from 
1752 to 1819.
}

\usage{data(Lynx)}
\format{
  A data frame with 68 observations on the following 2 variables.
  \describe{
    \item{\code{year}}{a numeric vector}
    \item{\code{pelts}}{a numeric vector}
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
data(Lynx)

xyplot(pelts ~ year, Lynx, type=c('p','l'))

\dontrun{
# Alternate form converting to Date class.
Year <- as.Date(paste("01jan", Lynx$date, sep = ""),
  "\%d\%b\%Y")
Lynx <- cbind(Lynx, Year)

require(ggplot2)
p <- ggplot(Lynx, aes(Year, no.pelts))
p + geom_line() + 
  geom_point(color = "red") +
  scale_y_continuous("Lynx fur returns") +
  opts(panel.grid.minor = theme_blank()) +
  opts(panel.grid.major = theme_line(size = 0.25, colour = "white"))}

}

\keyword{datasets}
