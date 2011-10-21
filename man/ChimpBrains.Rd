\name{ChimpBrains}    
\alias{ChimpBrains}
\docType{data}

\title{
Brodmann's Area 44 in Chimps
}

\description{
Asymmetry of Brodmann's area 44 in 20 chimpanzees.
}

\usage{data(ChimpBrains)}
\format{
  A data frame with 20 observations on the following 3 variables.
  \describe{
    \item{\code{name}}{name of chimp}
    \item{\code{sex}}{a factor with levels \code{F} and \code{M}}
    \item{\code{asymmetry}}{asymmetry score}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Cantalupo, C. and W.D. Hopkins. 2001. Asymmetric Broca's area in great apes. \emph{Nature} 414: 505.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(ChimpBrains)
xyplot(asymmetry ~ sex, ChimpBrains)
aggregate(asymmetry ~ sex, ChimpBrains, FUN = favstats)
if (require(Hmisc)) {
  summary(asymmetry ~ sex, ChimpBrains, fun = favstats)
}
}

\keyword{datasets}
