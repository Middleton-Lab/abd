\name{SockeyeFemaleBodyMass}
\alias{SockeyeFemaleBodyMass}
\docType{data}

\title{
Body masses of female Sockeye Salmon
}

\description{
Body Masses of 228 female Sockeye Salmon (\emph{Oncorhynchus nerka}; \url{http://www.nmfs.noaa.gov/pr/species/fish/sockeyesalmon.htm})
}

\usage{data(SockeyeFemaleBodyMass)}
\format{
  A data frame with 228 observations on the following variable.
  \describe{
    \item{\code{mass}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Hendry, A.P., O.K. Berg, and T.P. Quinn. 1999. Condition dependence and adaptation-by-time: Breeding date, life history, and energy allocation within a population of salmon. \emph{Oikos} 85: 499-514.
}

\references{
\url{http://www.jstor.org/stable/3546699}
}

\examples{
data(SockeyeFemaleBodyMass)

str(SockeyeFemaleBodyMass)
summary(SockeyeFemaleBodyMass)

dev.new(width = 9, height = 3)
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 3),
  xaxs = 'i',
  yaxs = 'i')
for (breaks in c(30, 10, 5)){
  hist(SockeyeFemaleBodyMass$mass, breaks = breaks,
  xlim = c(1, 4),
  col = 'red',
  ylab = 'Frequency',
  xlab = 'Body mass (kg)',
  main = '')
}
par(op)

}

\keyword{datasets}
