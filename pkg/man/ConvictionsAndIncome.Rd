\name{ConvictionsAndIncome}
\alias{ConvictionsAndIncome}
\docType{data}

\title{
ConvictionsAndIncome
}

\description{
Data reporting the number of individual with and without convictions per income level.
}

\usage{data(ConvictionsAndIncome)}
\format{
  A data frame with 6 observations on the following 3 variables.
  \describe{
    \item{\code{has.convictions}}{a factor with levels \code{no} and \code{yes}}
    \item{\code{income.level}}{a factor with levels \code{adequate}, \code{comfortable}, and \code{inadequate}}
    \item{\code{n}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Farrington, D.P. 1994. \emph{Cambridge Study in Delinquent Development} [Great Britain], 1961-1981. 2nd ICPSR ed. Inter-university Consortium for Political and Social Research, Ann Arbor, MI.
}

\references{
\url{http://webapp.icpsr.umich.edu/cocoon/NACJD-STUDY/08488.xml}
}

\examples{
data(ConvictionsAndIncome)
str(ConvictionsAndIncome)
ConvictionsAndIncome

# see https://stat.ethz.ch/pipermail/r-help/2009-January/185561.html
# for discussion of expand.dft(). Modified for ConvictionsAndIncome data.
expand.dft <- function(x, na.strings = 'NA', as.is = FALSE, dec = '.'){
  DF <- sapply(1:nrow(x), function(i) x[rep(i, each = x$n[i]), ], 
    simplify = FALSE)
  DF <- subset(do.call('rbind', DF), select = -n)
  for (i in 1:ncol(DF)){
    DF[[i]] <- type.convert(as.character(DF[[i]]),
      na.strings = na.strings, as.is = as.is, dec = dec)
  }
DF
} 

Conv.raw <- expand.dft(ConvictionsAndIncome)

xtabs(data = Conv.raw)

require(gmodels)
CrossTable(Conv.raw$has.convictions, Conv.raw$income.level,
  expected = FALSE,
  prop.r = FALSE,
  prop.c = FALSE,
  prop.chisq = FALSE, 
  prop.t = FALSE)
}

\keyword{datasets}
