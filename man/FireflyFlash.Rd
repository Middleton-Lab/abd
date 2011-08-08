\name{FireflyFlash}
\alias{FireflyFlash}
\docType{data}

\title{
Firefly Flash Duration
}

\description{
Flash duration (measured in milliseconds) of a sample of male fireflies (\emph{Photinus ignitus};
  \emph{n} = 35).
}

\usage{data(FireflyFlash)}
\format{
A data frame with 35 observations of one variable.
\describe{
\item{\code{flash}}{duration of flash (milliseconds)}
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Cratsley, C.K. and S.M. Lewis. 2003. Female preference
  for male courtship flashes in \emph{Photinus ignitus} fireflies.
  \emph{Behavioral Ecology} 14: 135-140. 
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(FireflyFlash)
str(FireflyFlash)
histogram(~flash, FireflyFlash)
}

\keyword{datasets}
