\docType{data}
\name{Guppies}
\alias{Guppies}
\title{Ornamentation and Attractiveness in Guppies}
\format{A data frame with 36 observations on the following 2 variables.
\describe{ \item{father.ornament}{a numeric vector}
\item{son.attract}{a numeric vector} }}
\source{
  \emph{inferred from} Brooks, R. 2000. Negative genetic
  correlation between male sexual attractiveness and
  survival. \emph{Nature} 406: 67-70.
}
\description{
  The father's ornamentation (composite score of color and
  brightness) and son's attractiveness (relative rates of
  visits by females) in male guppies (\emph{Poecilia
  reticulata}).
}
\examples{
data(Guppies)

str(Guppies)
xyplot(son.attract ~ father.ornament,
  Guppies,
  xlab = "Father's ornamentation",
  ylab = "Son's attractiveness"
  )
}
\references{
  \url{http://www.nature.com/nature/journal/v406/n6791/abs/406067a0.html}
}
\keyword{datasets}

