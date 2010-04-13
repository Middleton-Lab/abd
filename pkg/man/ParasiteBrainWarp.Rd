\name{ParasiteBrainWarp}
\alias{ParasiteBrainWarp}
\docType{data}

\title{
Frequencies of Fish Eaten by Trematode Infection Level
}

\description{
Frequencies of killifish (\emph{Fundulus parvipinnis}) eaten by birds depending on level of infection by the trematode \emph{Euhaplorchis californiensis}.
}

\usage{data(ParasiteBrainWarp)}
\format{
  A data frame with 6 observations on the following 3 variables.
  \describe{
    \item{\code{infection.status}}{a factor with levels: \code{high}, \code{light}, and \code{uninfected}}
    \item{\code{eaten}}{a factor with levels: \code{no} and \code{yes}}
    \item{\code{Frequency}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Lafferty, K.D. and A.K. Morris. 1996. Altered behavior of parasitized killifish increases susceptibility to predation by bird final hosts. \emph{Ecology} 77: 1390-1397.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(ParasiteBrainWarp)
ParasiteBrainWarp

ParasiteBrainWarp.long <- expand.dft(ParasiteBrainWarp, "Frequency")
table(ParasiteBrainWarp.long)
xtabs(~ eaten + infection.status, data = ParasiteBrainWarp.long)

require(gmodels)
CrossTable(ParasiteBrainWarp.long$eaten,
  ParasiteBrainWarp.long$infection.status,
  expected = TRUE,
  prop.r = FALSE,
  prop.c = FALSE,
  prop.chisq = TRUE, 
  prop.t = FALSE)
}


\keyword{datasets}
