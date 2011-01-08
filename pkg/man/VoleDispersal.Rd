\name{VoleDispersal}
\alias{VoleDispersal}
\docType{data}

\title{
Home Range Size in Field Voles
}

\description{
Home range size size in field voles (\emph{Microtus agrestis}).
}

\usage{data(VoleDispersal)}
\format{
  A data frame with 5 observations on the following 3 variables.
  \describe{
    \item{\code{homeranges}}{a numeric vector}
    \item{\code{count}}{a numeric vector}
    \item{\code{sex}}{a factor with levels \code{female} and \code{male}}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Sandell, M., J. Agrell, S. Erlinge, and J. Nelson. 1991. Adult philopatry and dispersal in the field vole \emph{Microtus agrestis}. \emph{Oecologia} 86: 153-158.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(VoleDispersal)
xtabs(count~sex+homeranges,VoleDispersal)
barchart( xtabs(count~sex+homeranges,VoleDispersal), auto.key=TRUE)
barchart(count~sex+homeranges,VoleDispersal)
barchart(count~sex,groups=homeranges,VoleDispersal)
barchart(count~sex,groups=homeranges,VoleDispersal,stack=TRUE)
}

\keyword{datasets}
