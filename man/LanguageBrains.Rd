\docType{data}
\name{LanguageBrains}
\alias{LanguageBrains}
\title{Brain Structure in Bilingual Humans}
\format{A data frame with 22 observations on the following 2 variables.
\describe{ \item{proficiency}{a numeric vector}
\item{greymatter}{a numeric vector} }}
\source{
  Mechelli, A., J.T. Crinion, U. Noppeney, J. O'Doherty, J.
  Ashburner, R.S. Frackowiak, and C.J. Price. 2004.
  Structural plasticity in the bilingual brain.
  \emph{Nature} 431: 757.
}
\description{
  Proficiency score (summary of reading, writing, and
  speech) in subjects' second language and density of gray
  matter in the left inferior parietal region.
}
\examples{
data(LanguageBrains)
str(LanguageBrains)
xyplot(proficiency ~ greymatter, LanguageBrains)
}
\keyword{datasets}

