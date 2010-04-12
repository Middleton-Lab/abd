\name{SpiderRunningAmputation}
\alias{SpiderRunningAmputation}
\docType{data}

\title{
Spider Running Speeds after Amputation
}

\description{
Data on speed before (\code{speed.before}) and speed after (\code{speed.after}) amputation of a pedipalp in the spider genus \emph{Tidarren}.
}

\usage{data(SpiderRunningAmputation)}
\format{
  A data frame with 32 observations on the following 2 variables.
  \describe{
    \item{\code{speed}}{a numeric vector}
    \item{\code{amputation.status}}{a factor with levels \code{after} and \code{before}}
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
data(SpiderRunningAmputation)

SpiderRunningAmputation
before <- subset(SpiderRunningAmputation, amputation.status == "before")
median(before$speed)

after <- subset(SpiderRunningAmputation, amputation.status == "after")
median(after$speed)

# Note that summary() returns quantiles.
(SRA.summary <- summary(before$speed))

# Interquantile range
SRA.summary[[5]] - SRA.summary[[2]]

# fivenum() produces quartiles, which matches the calculation on p. 69
SRA.5num <- fivenum(before$speed)
SRA.5num[[4]] - SRA.5num[[2]]

boxplot(speed ~ amputation.status, data = SpiderRunningAmputation,
  names = c("After amputation", "Before amputation"),
  ylab = "Running speed (cm/s)")

\dontrun{
# Using ggplot()
require(ggplot2)
p <- ggplot(SpiderRunningAmputation, aes(amputation.status, speed))
p + geom_boxplot() +
  scale_x_discrete("Amputation Status", 
    breaks = levels(SpiderRunningAmputation$amputation.status),
    labels = c("After amputation", "Before amputation")) +
  scale_y_continuous("Running speed (cm/s)")
}
}

\keyword{datasets}
