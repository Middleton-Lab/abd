\name{ConvictionsAndIncome}
\alias{ConvictionsAndIncome}
\docType{data}

\title{
ConvictionsAndIncome
}

\description{
ConvictionsAndIncome
}

\usage{data(ConvictionsAndIncome)}
\format{
  A data frame with 6 observations on the following 3 variables.
  \describe{
    \item{\code{has.convictions}}{a factor with levels \code{no} \code{yes}}
    \item{\code{income.level}}{a factor with levels \code{adequate} \code{comfortable} \code{inadequate}}
    \item{\code{n}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Farrington, D.P. 1994. Cambridge Study in Delinquent Development [Great Britain], 1961-1981. 2nd ICPSR ed. Inter-university Consortium for Political and Social Research, Ann Arbor, MI.
}

\references{
\url{http://webapp.icpsr.umich.edu/cocoon/NACJD-STUDY/08488.xml}
}

\examples{
data(ConvictionsAndIncome)
## maybe str(ConvictionsAndIncome) ; plot(ConvictionsAndIncome) ...
}

\keyword{datasets}
