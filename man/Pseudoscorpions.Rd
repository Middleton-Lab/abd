\name{Pseudoscorpions}
\alias{Pseudoscorpions}
\docType{data}

\title{
Multiple Mating in Pseudoscorpions
}

\description{
Successful numbers of broods (\code{Number.of.successful.broods}) in two groups of female pseudoscrpions (\emph{Cordylochernes scorpioides}), one mated to the same male twice and one to two different males.
}

\usage{data(Pseudoscorpions)}
\format{
  A data frame with 36 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{DM} \code{SM}}
    \item{\code{successful.broods}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Newcomer, S.D., J.A. Zeh, and D.W. Zeh. 1999. Genetic benefits enhance the reproductive success of polyandrous females. \emph{Proceedings of the National Academy of Sciences (USA)} 96: 10236-10241.
}

\references{
\url{http://www.pnas.org/content/96/18/10236.long}
}

\examples{
data(Pseudoscorpions)
str(Pseudoscorpions)
bwplot(successful.broods ~ treatment, Pseudoscorpions)
aggregate(successful.broods ~ treatment, Pseudoscorpions, favstats)
}

\keyword{datasets}
