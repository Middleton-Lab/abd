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
str(Aspirin.expanded)

# Plot 2X2 Contingency tables
plot( ~ treatment + cancer, data = Aspirin.expanded)
plot(table(Aspirin.expanded), main = "")

# Calculate odds
(Pr.asp <- 18496 / (18496 + 1438))
(Odds.asp <- Pr.asp / (1 - Pr.asp))
(Pr.no.asp <- 18515 / (18515 + 1427))
(Odds.no.asp <- Pr.no.asp / (1 - Pr.no.asp))
(Odds <- Odds.asp / Odds.no.asp)
ln.Odds <- log(Odds)

(SE.Odds <- sqrt(sum(1/Aspirin$count)))
Z <- 1.96
(CI.low <- ln.Odds - Z * SE.Odds)
(CI.high <- ln.Odds + Z * SE.Odds)

exp(CI.low)
exp(CI.high)

# Using odds.ratio()
# First reformat the data so that "No cancer" is in column 1
# and "Aspirin" is in row 2.
x <- matrix(c(18515, 18496, 1427, 1438), nrow = 2)
x
odds.ratio(x)
}

\keyword{datasets}
