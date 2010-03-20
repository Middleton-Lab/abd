\name{EelgrassGenotypes}
\alias{EelgrassGenotypes}
\docType{data}

\title{
Eelgrass Genotypes
}

\description{
Number of shoots (\code{shoots}) surviving in each of 32 experimental plots planted with 1, 3, or 6 different genotypes of eelgrass (\code{treatment.genotypes}).
}

\usage{data(EelgrassGenotypes)}
\format{
  A data frame with 32 observations on the following 2 variables.
  \describe{
    \item{\code{treatment.genotypes}}{a numeric vector of the number of genotypes planted in each plot}
    \item{\code{shoots}}{a numeric vector of the total number of shoots in each plot}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Reusch, T.B.H., A. Ehlers, A. Hämmerli, and B. Worm. 2005. Ecosystem recovery after climatic extremes enhanced by genotypic diversity. \emph{Proceedings of the National Academy of Sciences USA} 102: 2826-2831.
}

\references{
\url{http://www.pnas.org/content/102/8/2826.abstract}
}

\examples{
data(EelgrassGenotypes)
EelgrassGenotypes
EelgrassGenotypes$treatment.genotypes <-
  factor(EelgrassGenotypes$treatment.genotypes)
str(EelgrassGenotypes)
}

\keyword{datasets}
