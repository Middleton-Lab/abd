\docType{data}
\name{ProgesteroneExercise}
\alias{ProgesteroneExercise}
\title{Progesterone and Exercise}
\format{A data frame with 30 observations on the following 2 variables.
\describe{ \item{progesterone}{a numeric vector}
\item{ventilation}{a numeric vector} }}
\source{
  Brutsaert, T.D., H. Spielvogel, E. Caceres, M. Araoz,
  R.T. Chatterton, V.J. Vitzthum. 2002. Effect of menstrual
  cycle phase on exercise performance of high-altitude
  native women at 3600 m. \emph{Journal of Experimental
  Biology} 205: 233-239
}
\description{
  Progesterone levels and rates of ventilation during
  submaximal exercise in 30 women.
}
\examples{
data(ProgesteroneExercise)
str(ProgesteroneExercise)
xyplot(ventilation ~ progesterone, ProgesteroneExercise)
}
\references{
  \url{http://jeb.biologists.org/cgi/content/abstract/205/2/233}
}
\keyword{datasets}

