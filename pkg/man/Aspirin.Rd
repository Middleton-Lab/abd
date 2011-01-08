\name{Aspirin}
\alias{Aspirin}
\docType{data}

\title{
Effects of Aspirin on Cancer Rates
}

\description{
Frequency of cancer in 39,876 women taking and not taking aspirin.
}

\usage{data(Aspirin)}
\format{
  A data frame with 4 observations on the following 3 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{Aspirin} and \code{Placebo}}
    \item{\code{cancer}}{a factor with levels \code{no} and \code{yes}}
    \item{\code{count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Cook, N.R., I. Lee, J.M. Gaziano, D. Gordon, P.M. Ridker, J.E. Manson, C.H. Hennekens, and J.E. Buring. 2005. Low-dose aspirin in the primary prevention of cancer. \emph{Journal of the American Medical Association} 294: 47-55.
}

\references{
\url{http://jama.ama-assn.org/cgi/content/abstract/294/1/47}
}

\examples{
data(Aspirin)
Aspirin
Aspirin.expanded <- expand.dft(Aspirin, "count")
xtabs(~ cancer + treatment, Aspirin.expanded)
if (require(vcd)) {
	mosaic(~cancer + treatment, Aspirin.expanded)
}
}

\keyword{datasets}
