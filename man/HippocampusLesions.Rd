\docType{data}
\name{HippocampusLesions}
\alias{HippocampusLesions}
\title{Memory and the Hippocampus}
\format{A data frame with 57 observations on the following 2 variables.
\describe{ \item{lesion}{a numeric vector} \item{memory}{a
numeric vector} }}
\source{
  Broadbent, N.J., L.R. Squire, and R.E. Clark. 2004.
  Spatial memory, recognition memory, and the hippocampus.
  \emph{Proceedings of the National Academy of Sciences
  (USA)} 101: 14515-14520.
}
\description{
  Spatial memory score (\code{memory}) and percent lesion
  of the hippocampus (\code{lesion}).
}
\examples{
data(HippocampusLesions)
HippocampusLesions

xyplot(memory ~ lesion, data = HippocampusLesions,
  pch = 16, col = "red")

plot(memory ~ lesion, data = HippocampusLesions,
  pch = 16, col = "red")
}
\keyword{datasets}

