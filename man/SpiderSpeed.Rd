\name{SpiderSpeed}
\alias{SpiderSpeed}
\docType{data}

\title{
Spider Running Speeds after Amputation
}

\description{
Data on speed before and after amputation of a pedipalp in the 
spider genus \emph{Tidarren}.
}

\usage{data(SpiderSpeed)}
\format{
  A data frame with 32 observations on the following 2 variables.
  \describe{
    \item{\code{speed.before}}{speed (cm/s) before amputation }
    \item{\code{speed.after}}{speed (cm/s) after amputation }
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Ramos, M., D.J. Irschick, and T.E. Christenson. 2004. Overcoming an evolutionary conflict: Removal of a reproductive organ greatly increases locomotor performance. \emph{Proceedings of the National Academy of Sciences (USA)} 101: 4883-4887.
}

\references{
\url{http://www.pnas.org/content/101/14/4883}

\url{http://en.wikipedia.org/wiki/Pedipalp}

\url{http://en.wikipedia.org/wiki/Tidarren}
}

\examples{
data(SpiderSpeed)
xyplot(speed.after ~ speed.before, SpiderSpeed)
favstats(SpiderSpeed$speed.before)
favstats(SpiderSpeed$speed.after)
favstats(SpiderSpeed$speed.after - SpiderSpeed$speed.before)
}

\keyword{datasets}
