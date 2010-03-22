\name{SticklebackPlates}
\alias{SticklebackPlates}
\docType{data}

\title{
Number of Lateral Plates in Sticklebacks
}

\description{
Number of lateral plates (\code{no.plates}) in threespine sticklebacks (\emph{Gasterosteus aculeatus}) with three different \emph{Ectodysplasin} genotypes (\code{mm}, \code{Mm}, and \code{MM}).
}

\usage{data(SticklebackPlates)}
\format{
  A data frame with 344 observations on the following 2 variables.
  \describe{
    \item{\code{genotype}}{a factor with levels \code{mm}, \code{Mm}, and \code{MM}}
    \item{\code{no.plates}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Colosimo, P.F., C.L. Peichel, K. Nereng, B.K. Blackman, M.D. Shapiro, D. Schluter, and D.M. Kingsley. 2004. The genetic architecture of parallel armor plate reduction in threespine sticklebacks. \emph{PLOS Biology} 2: 635-641. \url{http://www.plosbiology.org/article/info:doi/10.1371/journal.pbio.0020109}}

\references{
Colosimo P.F., K.E. Hosemann, S. Balabhadra, G. Villarreal, M. Dickson, J. Grimwood, J Schmutz, R.M. Myers, D. Schluter, D.M. Kingsley. 2005. Widespread parallel evolution in sticklebacks by repeated fixation of ectodysplasin alleles. 307: 1928-33. \url{http://www.sciencemag.org/cgi/content/full/307/5717/1928}
}

\examples{
data(SticklebackPlates)

\dontrun{
op <- par(no.readonly = TRUE)
par(mfrow = c(3, 1),
  xaxs = "i",
  yaxs = "i")
for (i in c("mm", "Mm", "MM")){
  subset.by.genotype <- subset(SticklebackPlates, genotype == i)
  hist(subset.by.genotype$no.plates,
    breaks = 30,
    xlim = c(0, 70),
    ylim = c(0, 50),
    col = "red",
    ylab = "Frequency",
    xlab = "Number of Lateral Body Plates",
    main = paste(i))
}
par(op)

require(ggplot2)
p1 <- ggplot(SticklebackPlates, aes(no.plates))
p1 + geom_histogram(fill = "red", binwidth = 2) +
  facet_grid(genotype ~ .) +
  scale_x_continuous("Number of Lateral Body Plates") +
  scale_y_continuous("Frequency")
  
p2 <- ggplot(SticklebackPlates, aes(genotype, no.plates))
  p2 + geom_boxplot() +
  scale_x_discrete("Genotype") +
  scale_y_continuous("Number of Lateral Body Plates")
}

}

\keyword{datasets}
