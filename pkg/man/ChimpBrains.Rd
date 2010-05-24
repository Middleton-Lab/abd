\name{ChimpBrains}
\alias{ChimpBrains}
\docType{data}

\title{
Brodmann's Area 44 in Chimps
}

\description{
Asymmetry of Brodmann's area 44 in 20 chimpanzees.
}

\usage{data(ChimpBrains)}
\format{
  A data frame with 20 observations on the following 3 variables.
  \describe{
    \item{\code{Chimp.Name}}{a character vector}
    \item{\code{Sex}}{a factor with levels \code{F} and \code{M}}
    \item{\code{Rounded.AQ}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Cantalupo, C. and W.D. Hopkins. 2001. Asymmetric Broca's area in great apes. Nature 414: 505.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(ChimpBrains)
str(ChimpBrains)

# Bootstrap SE
set.seed(5)
n <- 10000
boot.median <- numeric(n)
for (i in 1:n){
  samp <- sample(ChimpBrains$Rounded.AQ, 20, replace = TRUE)
  boot.median[i] <- median(samp)
}
hist(boot.median, breaks = 50, col = "red")
mean(boot.median)
sd(boot.median) # Bootstrap standard error

\dontrun{
# Using boot()
require(boot)
boot.median <- function(x, i) median(x[i])
boot(ChimpBrains$Rounded.AQ, boot.median, 10000)}
}

\keyword{datasets}
