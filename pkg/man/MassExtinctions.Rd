\name{MassExtinctions}
\alias{MassExtinctions}
\docType{data}

\title{
Mass Extinction Frequency
}

\description{
The frequency of mass extinctions in the fossil record.
}

\usage{data(MassExtinctions)}
\format{
  A data frame with 21 observations on the following 2 variables.
  \describe{
    \item{\code{num.extinctions}}{a numeric vector}
    \item{\code{count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Raup, D.M. and J.J. Sepkoski, Jr. 1982. Mass extinctions in the marine fossil record. \emph{Science} 215: 1501-1503.
}

\references{
\url{http://www.sciencemag.org/cgi/content/abstract/sci;215/4539/1501}
}

\examples{
data(MassExtinctions)
MassExtinctions

if(0){
# Calculate weighted mean
# with expand.dft()
n.extinctions <- expand.dft(MassExtinctions,
  "count")$Number.of.extinctions
wt.mean <- mean(n.extinctions)

# With weighted.mean()
(wt.mean <- weighted.mean(MassExtinctions$Number.of.extinctions,
  MassExtinctions$Frequency))

hist(n.extinctions,
  ylim = c(0, 30),
  xlab = "Number of Extinctions",
  main = "Frequency of Mass Extinctions")

(Pr.3 <- (exp(-wt.mean) * wt.mean^3) / factorial(3))
76 * Pr.3

# With dpois()
76 * dpois(3, wt.mean)

# Calculate expected
expected <- (exp(-wt.mean) * wt.mean^c(0:20) /
  factorial(c(0:20))) * 76
# Alternately with dpois()
expected.dpois <- dpois((0:20), lambda = wt.mean) * 76

# Compare the two
data.frame(expected, expected.dpois)

# Collapse some rows into a single expected value
expected2 <- c(sum(expected[1:2]), expected[3:8], 
  sum(expected[9:21]))
expected2

MassExtinctions2 <- rbind(MassExtinctions[-c(1, 9:21), ], c(8, 9))
MassExtinctions2

chisq <- sum((MassExtinctions2$Frequency - expected2)^2 / 
  expected2)
chisq
1 - pchisq(chisq, df = 6)

# Alternate using chisq.test(). Note that df = 7 here, because 
# chisq.test() doesn't know that mu was estimated from the data.
chisq.test(MassExtinctions2$Frequency, p = expected2,
  rescale.p = TRUE)

\dontrun{
# Second alternate using goodfit() from vcd package
require(vcd)
extinctions.fit <- goodfit(MassExtinctions2$Frequency)
summary(extinctions.fit)}

# Variance
var(n.extinctions)
}
}

\keyword{datasets}
