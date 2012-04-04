\docType{data}
\name{SticklebackPlates}
\alias{SticklebackPlates}
\title{Number of Lateral Plates in Sticklebacks}
\format{A data frame with 344 observations on the following 2 variables.
\describe{ \item{genotype}{a factor with levels \code{mm},
\code{Mm}, and \code{MM}} \item{plates}{number of plates} }}
\source{
  Colosimo, P.F., C.L. Peichel, K. Nereng, B.K. Blackman,
  M.D. Shapiro, D. Schluter, and D.M. Kingsley. 2004. The
  genetic architecture of parallel armor plate reduction in
  threespine sticklebacks. \emph{PLoS Biology} 2: 635-641.
  \url{http://www.plosbiology.org/article/info:doi/10.1371/journal.pbio.0020109}
}
\description{
  Number of lateral plates (\code{plates}) in threespine
  sticklebacks (\emph{Gasterosteus aculeatus}) with three
  different \emph{Ectodysplasin} genotypes (\code{mm},
  \code{Mm}, and \code{MM}).
}
\examples{
data(SticklebackPlates)

aggregate(plates ~ genotype, SticklebackPlates, FUN = favstats)

histogram( ~ plates | genotype, SticklebackPlates,
  layout = c(1,3),
  n = 15,
  xlab = "Number of Lateral Body Plates"
  )

densityplot( ~ plates | genotype, SticklebackPlates,
  xlab = "Number of Lateral Body Plates",
  layout = c(1,3)
  )
}
\references{
  Colosimo P.F., K.E. Hosemann, S. Balabhadra, G.
  Villarreal, M. Dickson, J. Grimwood, J Schmutz, R.M.
  Myers, D. Schluter, D.M. Kingsley. 2005. Widespread
  parallel evolution in sticklebacks by repeated fixation
  of ectodysplasin alleles. \emph{Science }307: 1928-33.
  \url{http://www.sciencemag.org/cgi/content/full/307/5717/1928}
}
\keyword{datasets}

