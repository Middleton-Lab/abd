\name{SockeyeFemaleBodyMass}
\alias{SockeyeFemaleBodyMass}
\docType{data}

\title{
Body Masses of Female Sockeye Salmon
}

\description{
Body Masses of Female Sockeye Salmon
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
Hendry, A.P., O.K. Berg, and T.P. Quinn. 1999. Condition Dependence and Adaptation-by-Time: Breeding Date, Life History, and Energy Allocation within a Population of Salmon. \emph{Oikos} 85: 499-514.
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
