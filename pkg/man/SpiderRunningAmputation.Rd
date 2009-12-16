\name{SpiderRunningAmputation}
\alias{SpiderRunningAmputation}
\docType{data}

\title{
SpiderRunningAmputation
}

\description{
Data on speed before (\code{speed.before}) and speed after (\code{speed.after}) amputation of a pedipalp (\url{http://en.wikipedia.org/wiki/Pedipalp}) in the spider genus \emph{Tidarren} (\url{http://en.wikipedia.org/wiki/Tidarren}).
}

\usage{data(SpiderRunningAmputation)}
\format{
  A data frame with 16 observations on the following 2 variables.
  \describe{
    \item{\code{speed.before}}{a numeric vector}
    \item{\code{speed.after}}{a numeric vector}
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
}

\examples{
data(SpiderRunningAmputation)

SpiderRunningAmputation
median(SpiderRunningAmputation$speed.before)

# Note that the values for the 1st and 3rd quartiles reported
# differ from those given on p. 68. See footnote 5.
(SRA.summary <- summary(SpiderRunningAmputation$speed.before))

# Interquartile range
SRA.summary[[5]] - SRA.summary[[2]]

# Reformat data to 'long' format for boxplots
SRA.long <- stack(SpiderRunningAmputation)
names(SRA.long) <- c('speed', 'amputation.status')
boxplot(speed ~ amputation.status, data = SRA.long,
  names = c('After amputation', 'Before amputation'),
  ylab = 'Running speed (cm/s)')

\dontrun{
# Using ggplot()
require(ggplot2)
p <- ggplot(SRA.long, aes(amputation.status, speed))
p + geom_boxplot() +
  scale_x_discrete('', labels = c('After amputation', 'Before amputation')) +
  scale_y_continuous('Running speed (cm/s)')}

}

\keyword{datasets}
