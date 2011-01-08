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
    \item{\code{genotype}}{a factor with levels \code{F1} and \code{F2}}
    \item{\code{preference}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Haeslery, M.P. and O. Seehausen. 2005. Inheritance of female mating preference
in a sympatric sibling species pair of Lake Victoria cichlids: implications for
speciation. \emph{Proceedings of the Royal Society of London, Series B,
Biological Sciences} 272: 237-245.
}

\references{
\url{http://rspb.royalsocietypublishing.org/content/272/1560/237.full.pdf}
}

\examples{
data(Cichlids)
str(Cichlids)

if (require(Hmisc)) {
	summary(preference ~ genotype, Cichlids, fun=favstats)
} else {
	aggregate(preference ~ genotype, Cichlids, FUN=favstats)
}

if (require(plyr)) {
ddply(Cichlids, .(genotype),
  function(df)c(mean = mean(df$preference),
                standard.deviation = sd(df$preference),
                n = length(df$preference)))
}
}

\keyword{datasets}
