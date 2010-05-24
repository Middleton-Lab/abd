\name{BrookTrout}
\alias{BrookTrout}
\docType{data}
\title{
Salmon Survival in the Presence of Brook Trout
}

\description{
Total numbers of salmon released (\code{salmon.released}) and surviving (\code{salmon.surviving}) in 12 streams, 6 with brook trout \code{present} and 6 with brook trout \code{absent}. The proportion of salmon surviving (\code{proportion.surviving}) is given for each stream.
}

\usage{data(BrookTrout)}

\format{
  A data frame with 12 observations on the following 4 variables.
  \describe{
    \item{\code{brook.trout}}{a factor with levels \code{absent} and \code{present} indicating whether brook trout are absent or present in the stream}
    \item{\code{salmon.released}}{a numeric vector of the total number of salmon released}
    \item{\code{salmon.surviving}}{a numeric vector of the number of salmon surviving}
    \item{\code{proportion.surviving}}{a numeric vector of the proportion of salmon surviving}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Levin, P.S., S. Achord, B.E. Fiest, and R.W. Zabel. 2002. Non-indigenous brook trout and the demise of Pacific salmon: a forgotten threat? \emph{Proceedings of the Royal Society of London, Series B, Biological Sciences} 269: 1663-1670.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(BrookTrout)
str(BrookTrout)
BrookTrout

# Aggregate the data using ddply()
require(plyr)
salmon.aggregate <- ddply(BrookTrout, .(brook.trout),
  function(x)c(sum(x$salmon.released - x$salmon.surviving),
    sum(x$salmon.surviving)))
names(salmon.aggregate)[c(2,3)] <- c("Survived", "Died")
salmon.aggregate

# Boxplot
boxplot(proportion.surviving ~ brook.trout, data = BrookTrout,
  ylab = "Proportion Surviving",
  names = c("Trout Absent", "Trout Present"))

\dontrun{
# Dotplot
require(lattice)
dotplot(proportion.surviving ~ brook.trout, data = BrookTrout)}

# Aggregate again, calculating mean, standard deviation, and n
require(plyr)
salmon.aggregate2 <- ddply(BrookTrout, .(brook.trout),
  function(x)c(mean(x$proportion.surviving),
               sd(x$proportion.surviving),
               length(x$proportion.surviving)))
names(salmon.aggregate2) <- c("Group", "Sample Mean",
  "Sample Standard Deviation", "Sample Size")
salmon.aggregate2

# Use Welch's t-test, because the variances are not equal
t.test(proportion.surviving ~ brook.trout, data = BrookTrout,
  var.equal = FALSE)

# Comparing variances
# Bartlett Test
bartlett.test(proportion.surviving ~ brook.trout, data = BrookTrout)

# Exploring differences in variance
set.seed(2)
x1 <- rnorm(40, sd = 1)
x2 <- rnorm(40, sd = 1)

x12 <- c(x1, x2)
A <- factor(rep(1:2, each = 40))

plot(density(x12), type = "n", ylim = c(0, 0.5))
lines(density(x1), col = "blue")
lines(density(x2), col = "red")

bartlett.test(x12 ~ A)

# Same Mean; Different sd
x2 <- rnorm(40, sd = 2)
x12 <- c(x1, x2)
dev.new()
plot(density(x12), type = "n", ylim = c(0, 0.5))
lines(density(x1), col = "blue")
lines(density(x2), col = "red")
bartlett.test(x12 ~ A)
}

\keyword{datasets}
