\docType{data}
\name{ChimpBrains}
\alias{ChimpBrains}
\title{Brodmann's Area 44 in Chimps}
\format{A data frame with 20 observations on the following 3 variables.
\describe{ \item{name}{name of chimp} \item{sex}{a factor
with levels \code{F} and \code{M}} \item{asymmetry}{asymmetry score}
}}
\source{
  Cantalupo, C. and W.D. Hopkins. 2001. Asymmetric Broca's
  area in great apes. \emph{Nature} 414: 505.
}
\description{
  Asymmetry of Brodmann's area 44 in 20 chimpanzees.
}
\examples{
data(ChimpBrains)
xyplot(asymmetry ~ sex, ChimpBrains)
aggregate(asymmetry ~ sex, ChimpBrains, FUN = favstats)
summary(asymmetry ~ sex, ChimpBrains, fun = favstats)
}
\keyword{datasets}

