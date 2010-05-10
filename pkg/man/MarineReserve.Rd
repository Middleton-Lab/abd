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
A numeric vector with 32 observations.
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Halpern, B. S. 2003. The impact of marine reserves: do reserves work and does reserve size matter? \emph{Ecological Applications} 13: S117-S137.
}

\references{
\url{http://www.nceas.ucsb.edu/~halpern/pdf/Halpern_EA_2003.pdf}
}

\examples{
data(MarineReserve)
str(MarineReserve)

hist(MarineReserve)

# Normal quantile plot; Note that the default is datax = FALSE
qqnorm(MarineReserve, datax = TRUE)
qqline(MarineReserve, datax = TRUE)

# Natural log transformation
log.biomass <- log(MarineReserve)
hist(log.biomass)
(mean(log.biomass))
(sd(log.biomass))

t.test(log.biomass, mu = 0, var.equal = TRUE)

# Confidence intervals
(cis <- ci(log.biomass))

# Back transform
exp(cis$lower)
exp(cis$upper)
}

\keyword{datasets}
