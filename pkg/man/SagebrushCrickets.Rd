\name{SagebrushCrickets}
\alias{SagebrushCrickets}
\docType{data}

\title{
Sagebrush Cricket Mating Times
}

\description{
Time to mating (\code{Time.to.mating}) in fed and unfed (\code{Feeding}) sagebrush crickets (\emph{Cyphoderris strepitans}).
}

\usage{data(SagebrushCrickets)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{\code{Feeding}}{a factor with levels \code{fed} \code{starved}}
    \item{\code{Time.to.mating}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Chadwick Johnson, J., T.M. Ivy, and S.K. Sakaluk. 1999. Female remating propensity contingent on sexual cannibalism in sagebrush crickets, \emph{Cyphoderris strepitans}: a mechanism of cryptic female choice. \emph{Behavioral Ecology} 10: 227-233.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(SagebrushCrickets)
SagebrushCrickets
str(SagebrushCrickets)

# Subset and extract the Time.to.mating data
starved <- subset(SagebrushCrickets,
  Feeding == "starved")$Time.to.mating
fed <- subset(SagebrushCrickets, Feeding == "fed")$Time.to.mating

dev.new()
par(mfrow = c(2, 1))
hist(starved, xlim = c(0, 100))
hist(fed, xlim = c(0, 100))

# Sort the SagebrushCrickets data.frame
sorted <- SagebrushCrickets[order(SagebrushCrickets$Time.to.mating), ]

# Add a rank column
sorted$rank <- 1:24
sorted

# Extract n
(n.fed <- length(fed))
(n.starved <- length(starved))

# Calculate rank sum
(sum.fed <- sum(sorted$rank[sorted$Feeding == "fed"]))
(sum.starved <- sum(sorted$rank[sorted$Feeding == "starved"]))

# Calculate U for each group
(u.starved <- n.starved * n.fed + 
  (n.starved * (n.starved + 1) / 2) - sum.starved)
(u.fed <- n.fed * n.starved - u.starved)

# Choose the larger U
(u <- max(c(u.starved, u.fed)))

# Critical value for p = 0.05, with n1 = 11 and n2 = 13
qwilcox(1-(0.05/2), 11, 13)

# Alternately with wilcox.test()
wilcox.test(Time.to.mating ~ Feeding, data = SagebrushCrickets)
}

\keyword{datasets}
