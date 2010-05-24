\name{Pseudoscorpions}
\alias{Pseudoscorpions}
\docType{data}

\title{
Multiple Mating in Pseudoscorpions
}

\description{
Successful numbers of broods (\code{Number.of.successful.broods}) in two groups of female pseudoscrpions (\emph{Cordylochernes scorpioides}), one mated to the same male twice and one to two different males.
}

\usage{data(Pseudoscorpions)}
\format{
  A data frame with 36 observations on the following 2 variables.
  \describe{
    \item{\code{Mating.treatment}}{a factor with levels \code{DM} \code{SM}}
    \item{\code{Number.of.successful.broods}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Newcomer, S.D., J.A. Zeh, and D.W. Zeh. 1999. Genetic benefits enhance the reproductive success of polyandrous females. \emph{Proceedings of the National Academy of Sciences (USA)} 96: 10236-10241.
}

\references{
\url{http://www.pnas.org/content/96/18/10236.long}
}

\examples{
data(Pseudoscorpions)
str(Pseudoscorpions)

# Shorten names
PS <- Pseudoscorpions
names(PS) <- c("tx", "broods")

obs.diff <- mean(PS$broods[PS$tx == "SM"]) - 
  mean(PS$broods[PS$tx == "DM"])
obs.diff

set.seed(6)
n <- 10000
diffs <- numeric(n)
for (i in 1:n){
  rand.broods <- sample(PS$broods)
  diffs[i] <- mean(rand.broods[PS$tx == "SM"]) - 
    mean(rand.broods[PS$tx == "DM"])
}

sum(diffs < obs.diff)/n

hist(diffs, breaks = 50)
abline(v = obs.diff, col = "red")

# As a two-sample t-test
# (much slower)
obs.t <- t.test(broods ~ tx, data = PS)$statistic
set.seed(6)
n <- 10000
ts <- numeric(n)
for (i in 1:n){
  rand.broods <- sample(PS$broods)
  ts[i] <- t.test(rand.broods ~ tx, data = PS)$statistic
}

sum(ts > obs.t)/n

hist(ts, breaks = 50)
abline(v = obs.t, col = "red")
}

\keyword{datasets}
