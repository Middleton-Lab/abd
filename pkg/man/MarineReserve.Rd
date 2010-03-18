\name{MarineReserve}
\alias{MarineReserve}
\docType{data}

\title{
Marine Reserve Biomass
}

\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}

\usage{data(MarineReserve)}
\format{
  The format is:
 num [1:32] 1.34 1.96 2.49 1.27 1.19 1.15 1.29 1.05 1.1 1.21 ...
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

# Compare with a normal distribution of the same n
#  with the same mean
set.seed(3)
y <- rnorm(length(MarineReserve), mean = mean(MarineReserve))
qqnorm(y, datax = TRUE)
qqline(y, datax = TRUE)

# Natural log transformation
log.biomass <- log(MarineReserve)
hist(log.biomass)
(mean(log.biomass))
(sd(log.biomass))

t.test(log.biomass, mu = 0, var.equal = TRUE)

# Confidence intervals
require(gmodels)
(cis <- ci(log.biomass))

# Back transform
exp(cis[2])
exp(cis[3])
}

\keyword{datasets}
