\docType{data}
\name{Hemoglobin}
\alias{Hemoglobin}
\title{Hemoglobin Levels in High Altitude Populations}
\format{A data frame with 40 observations on the following 3 variables.
\describe{ \item{hemoglobin}{a numeric vector}
\item{group}{a factor with levels: \code{Andes}, \code{Ethiopia},
\code{Tibet}, and \code{USA}} \item{relative.frequency}{a numeric
vector} }}
\source{
  \emph{inferred from} Beall, C.M., M.J. Decker, G.M.
  Bittenham, I. Kushner, A. Gebremedhin, K.P. Strohl. 2002.
  An Ethiopian pattern of human adaptation to high-altitude
  hypoxia. \emph{Proceeding of the National Academy of
  Sciences (USA)} 99(26): 17215-17218.
}
\description{
  Relative rates of hemoglobin concentration in four
  populations of humans living at different altitudes.
}
\examples{
data(Hemoglobin)
str(Hemoglobin)

xyplot(relative.frequency ~ hemoglobin | group, Hemoglobin,
  type ='h', lwd=4, layout=c(1,4))
}
\references{
  \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC139295/}
}
\keyword{datasets}

