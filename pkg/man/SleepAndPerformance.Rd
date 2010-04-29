\name{SleepAndPerformance}
\alias{SleepAndPerformance}
\docType{data}

\title{
Sleep and Learning
}

\description{
The increase in "slow-wave" sleep and improvements in spatial learning tasks in 10 humans.
}

\usage{data(SleepAndPerformance)}
\format{
  A data frame with 10 observations on the following 2 variables.
  \describe{
    \item{\code{sleep}}{a numeric vector}
    \item{\code{improvement}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Huber, R., M.F. Ghilardi, M. Massimini, and G. Tononi. 2004. Local sleep and learning. \emph{Nature} 430: 78-81.
}

\references{
\url{http://www.ncbi.nlm.nih.gov/pubmed/15184907}

\url{http://www.nature.com/nature/journal/v430/n6995/full/nature02663.html}
}

\examples{
data(SleepAndPerformance)
str(SleepAndPerformance)
SleepAndPerformance
}

\keyword{datasets}
