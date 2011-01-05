\name{BirdSexRatio}
\alias{BirdSexRatio}
\docType{data}
\title{
Sex Ratios in Birds
}

\description{
Correlation coefficient of sex ratio in bird offspring.
}

\usage{data(BirdSexRatio)}
\format{
A data frame with 15 observations of one variable
\describe{
\item{\code{corr.coeff}}{correlation coefficient of sex ratio in bird offspring}
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
West, S.A. and B.C. Sheldon. 2002. Constraints in the evolution of sex ratio adjustment. \emph{Science} 295: 1695-1688.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(BirdSexRatio)
histogram(~corr.coeff, BirdSexRatio, n = 10,
  xlab = "Correlation Coefficient")
}

\keyword{datasets}
