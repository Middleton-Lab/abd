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

hist(MarineReserve$biomass.ratio)

# Normal quantile plot; Note that the default is datax = FALSE
qqmath(~biomass.ratio, MarineReserve)
qqnorm(MarineReserve$biomass.ratio, datax = TRUE)
qqline(MarineReserve$biomass.ratio, datax = TRUE)

# Natural log transformation
log.biomass <- log(MarineReserve$biomass)
hist(log.biomass)
(mean(log.biomass))
(sd(log.biomass))

t.test(log.biomass, mu = 0, var.equal = TRUE)

# Confidence intervals
cis <- interval( t.test(log.biomass, mu = 0, var.equal = TRUE) )

# Back transform
exp(cis)
}

\keyword{datasets}
