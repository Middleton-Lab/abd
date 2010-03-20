\name{DungBeetleCondition}
\alias{DungBeetleCondition}
\docType{data}

\title{
Heritability of Body Condition in Dung Beetles
}

\description{
Body condition (\code{offspring.condition}) in 36 dung beetles (\emph{Onthophagus taurus}) from 12 \code{male}s each mated to 3 different virgin females.
}

\usage{data(DungBeetleCondition)}
\format{
  A data frame with 36 observations on the following 2 variables.
  \describe{
    \item{\code{male}}{a numeric vector}
    \item{\code{offspring.condition}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Kotiaho, J.S., L.W. Simmons, and J.L. Tomkins. 2001. Towards a resolution of the lek paradox. \emph{Nature} 410: 684-686.
}

\references{
\url{http://en.wikipedia.org/wiki/Dung_beetle}

\url{http://www.nature.com/nature/journal/v410/n6829/abs/410684a0.html}
}

\examples{
data(DungBeetleCondition)
str(DungBeetleCondition)

DungBeetleCondition$male <- factor(DungBeetleCondition$male)
str(DungBeetleCondition)
}

\keyword{datasets}
