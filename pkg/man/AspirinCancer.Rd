\name{AspirinCancer}
\alias{AspirinCancer}
\docType{data}

\title{
Effects of Aspirin on Cancer Rates
}

\description{
Frequency of cancer in 39,876 women taking and not taking aspirin.
}

\usage{data(AspirinCancer)}
\format{
  A data frame with 4 observations on the following 3 variables.
  \describe{
    \item{\code{Aspirin.treatment}}{a factor with levels \code{Aspirin} \code{Placebo}}
    \item{\code{Cancer}}{a factor with levels \code{no} \code{yes}}
    \item{\code{Frequency}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Cook, N.R., I. Lee, J.M. Gaziano, D. Gordon, P.M. Ridker, J.E. Manson, C.H. Hennekens, and J.E. Buring. 2005. Low-dose aspirin in the primary prevention of cancer. Journal of the American Medical Association 294: 47-55.
}

\references{
\url{http://jama.ama-assn.org/cgi/content/abstract/294/1/47}
}

\examples{
data(AspirinCancer)
AspirinCancer

AspirinCancer.expanded <- expand.dft(AspirinCancer)
str(AspirinCancer.expanded)

# Plot 2X2 Contingency tables
plot( ~ Aspirin.treatment + Cancer, data = AspirinCancer.expanded)
plot(table(AspirinCancer.expanded), main = "")

# Calculate odds
(Pr.asp <- 18496 / (18496 + 1438))
(Odds.asp <- Pr.asp / (1 - Pr.asp))
(Pr.no.asp <- 18515 / (18515 + 1427))
(Odds.no.asp <- Pr.no.asp / (1 - Pr.no.asp))
Odds.asp / Odds.no.asp
}

\keyword{datasets}
