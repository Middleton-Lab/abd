\name{SockeyeFemaleBodyMass}
\alias{SockeyeFemaleBodyMass}
\docType{data}

\title{
Body Masses of Female Sockeye Salmon
}

\description{
Body Masses of 228 female Sockeye Salmon (\emph{Oncorhynchus nerka}; \url{http://www.nmfs.noaa.gov/pr/species/fish/sockeyesalmon.htm})
}

\usage{data(SockeyeFemaleBodyMass)}
\format{
  The format is:
 num [1:228] 3.09 2.91 3.06 2.69 2.88 2.98 1.61 2.16 1.56 1.76 ...
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
  xaxs = "i",
  yaxs = "i")
for (breaks in c(30, 10, 5)){
  hist(SockeyeFemaleBodyMass, breaks = breaks,
    xlim = c(1, 4),
    col = "red",
    ylab = "Frequency",
    xlab = "Body mass (kg)",
    main = "")
}

par(op)
}

\keyword{datasets}
