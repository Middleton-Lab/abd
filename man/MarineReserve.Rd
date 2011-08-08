\name{MarineReserve}
\alias{MarineReserve}
\docType{data}

\title{
Marine Reserve Biomass
}

\description{
Relative biomass in 32 marine reserves.
}

\usage{data(MarineReserve)}
\format{
A data frame with
32 observations of one variable.  
\describe{
\item{\code{biomass.ratio}}{}
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Halpern, B.S. 2003. The impact of marine reserves: do reserves work and does reserve size matter? \emph{Ecological Applications} 13: S117-S137.
}

\references{
\url{http://www.nceas.ucsb.edu/~halpern/pdf/Halpern_EA_2003.pdf}
}

\examples{
data(MarineReserve)
str(MarineReserve)
histogram(~biomass.ratio, MarineReserve)
}

\keyword{datasets}
