\name{Iguanas}
\alias{Iguanas}
\docType{data}

\title{
Iguana Body Length Changes
}

\description{
Body size change in 64 Galápagos marine iguanas (\emph{Amblyrhynchus cristatus}) that survived the 1992-1993 El Niño event.
}

\usage{data(Iguanas)}

\format{
A numeric vector with 64 observations.
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Wikelski, M. and C. Thom. 2000. Marine iguanas shrink to survive El Niño. \emph{Nature} 403: 37-38.
}

\references{
\url{http://en.wikipedia.org/wiki/Marine_iguana}
}

\examples{
data(Iguanas)
str(Iguanas)
hist(Iguanas, breaks = 10)
}

\keyword{datasets}
