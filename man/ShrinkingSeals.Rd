\docType{data}
\name{ShrinkingSeals}
\alias{ShrinkingSeals}
\title{Seal Body Lengths and Age}
\format{A data frame with 9,665 observations on the following 2 variables.
\describe{ \item{age}{age (days)} \item{length}{body length
(cm) } }}
\source{
  Trites, A.W. and M.A. Bigg. 1996. Physical growth of
  northern fur seals: seasonal fluctuations and migratory
  influences. \emph{Journal of Zoology (London)} 238:
  459-482.
}
\description{
  Body length (cm) and age (days) for 9,665 female Northern
  fur seals (\emph{Callorhinus ursinus}).
}
\examples{
data(ShrinkingSeals)
str(ShrinkingSeals)

plot(ShrinkingSeals, pch = 16, cex = 0.5)
xyplot(length ~ age, ShrinkingSeals, pch=16, alpha=.65, cex=.6)
}
\keyword{datasets}

