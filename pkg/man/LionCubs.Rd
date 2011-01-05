\name{LionCubs}
\alias{LionCubs}
\docType{data}

\title{
Time to Reproduction in Female Lions
}

\description{
Time to reproduction (\code{Days}) based on whether death of previous cubs
was due to infanticide (\code{New}) or accidental (\code{Same}).
}

\usage{data(LionCubs)}
\format{
  A data frame with 14 observations on the following 2 variables.
  \describe{
    \item{\code{cause.of.death}}{a factor with \code{accident} and \code{infanticide}}
    \item{\code{days.to.next.cub}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Packer, C. and A.E. Pusey. 1983. Adaptations of female lions to infanticide by incoming males. \emph{The American Naturalist} 121: 716-728.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(LionCubs)
xyplot(days.to.next.cub ~ cause.of.death, LionCubs)
}

\keyword{datasets}
