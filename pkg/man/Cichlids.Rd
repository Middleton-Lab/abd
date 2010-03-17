\name{Cichlids}
\alias{Cichlids}
\docType{data}

\title{
Cichlid Mating Preference
}

\description{
Preference index in F1 and F2 crosses of two species of cichlids from Lake Victoria, \emph{Pundamilia pundamilia} and \emph{P. nyererei}.
}

\usage{data(Cichlids)}
\format{
  A data frame with 53 observations on the following 2 variables.
  \describe{
    \item{\code{Genotype}}{a factor with levels \code{F1} \code{F2}}
    \item{\code{preference}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Haeslery, M.P. and O. Seehausen. 2005. \emph{Proceedings of the Royal Society of London, Series B, Biological Sciences} 272: 237-245.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(Cichlids)
str(Cichlids)

require(plyr)
ddply(Cichlids, .(Genotype),
  function(df)c(mean = mean(df$preference),
                standard.deviation = sd(df$preference),
                n = length(df$preference)))
}

\keyword{datasets}
