\name{TsetseLearning}
\alias{TsetseLearning}
\docType{data}

\title{
Dietary Learning in Tsetse Flies
}

\description{
Dietary conditioning \code{treatment} and subsequent proportion of tsetse flies (\emph{Glossina palpalis}) feeding on cow blood in each of 13 cohorts.
}

\usage{data(TsetseLearning)}
\format{
  A data frame with 13 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{cow} and \code{lizard}}
    \item{\code{proportion.cow}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Bouyer, J., M. Pruvot, Z. Bengaly, P.M. Guerin, and R. Lancelot. 2007. Learning influences host choice in tsetse. \emph{Biology Letters} 3: 113-116.
}

\references{
\url{http://rsbl.royalsocietypublishing.org/content/3/2/113.full}
}

\examples{
data(TsetseLearning)
TsetseLearning
str(TsetseLearning)
}

\keyword{datasets}
