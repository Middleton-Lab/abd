\name{BrookTrout}
\alias{BrookTrout}
\alias{BrookTrout2}
\docType{data}
\title{
Salmon Survival in the Presence of Brook Trout
}

\description{
Total numbers of salmon released (\code{salmon.released}) and surviving
(\code{salmon.surviving}) in 12 streams, 6 with brook trout \code{present} and
6 with brook trout \code{absent}. The proportion of salmon surviving
(\code{proportion.surviving}) is given for each stream.
}

\usage{
data(BrookTrout)

data(BrookTrout2)
}

\format{
  \code{BrookTrout} is a data frame with 12 observations on the following 
  4 variables.  \code{BrookTrout2} is a different summary of the same
  study and gives survival rates for chinook in different years.
  \describe{
    \item{\code{trout}}{a factor with levels \code{absent} and \code{present} indicating whether brook trout are absent or present in the stream}
    \item{\code{salmon.released}}{a numeric vector of the total number of salmon released}
    \item{\code{salmon.surviving}}{a numeric vector of the number of salmon surviving}
    \item{\code{proportion.surviving}}{a numeric vector of the proportion of salmon surviving}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Levin, P.S., S. Achord, B.E. Fiest, and R.W. Zabel. 2002. Non-indigenous brook trout and the demise of Pacific salmon: a forgotten threat? \emph{Proceedings of the Royal Society of London, Series B, Biological Sciences} 269: 1663-1670.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(BrookTrout)
data(BrookTrout2)
str(BrookTrout)
str(BrookTrout2)
bwplot(proportion.surviving ~ trout, BrookTrout)
aggregate( proportion.surviving ~ trout, BrookTrout, FUN = favstats)
if (require(Hmisc)) {
	summary( proportion.surviving ~ trout, BrookTrout, fun = favstats)
}
}

\keyword{datasets}
