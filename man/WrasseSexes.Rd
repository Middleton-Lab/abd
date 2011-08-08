\name{WrasseSexes}
\alias{WrasseSexes}
\docType{data}

\title{
Distribution of Wrasses
}

\description{
Number and sex of adult wrasses in a section of the Great Barrier Reef.
}

\usage{data(WrasseSexes)}
\format{
  A data frame with 3 observations on the following 3 variables.
  \describe{
    \item{\code{males}}{a numeric vector}
    \item{\code{females}}{a numeric vector}
    \item{\code{count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(WrasseSexes)
xtabs(count ~ males + females, WrasseSexes)
}

\keyword{datasets}
