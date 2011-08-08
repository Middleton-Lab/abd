\name{WillsPresidents}
\alias{WillsPresidents}
\docType{data}

\title{
Presidential "Wills"
}

\description{
Number of times a presidential candidate said "will," "shall," or "going to" in presidential debates from 1960-2004 (years incomplete).
}

\usage{data(WillsPresidents)}
\format{
  A data frame with 16 observations on the following 3 variables.
  \describe{
    \item{\code{candidate}}{a character vector with the candidate's name}
    \item{\code{winner}}{a factor with levels \code{n} \code{y} indicating whether the candidate won the election \code{y} or not.}
    \item{\code{wills}}{a numeric vector}
    \item{\code{loser.wills}}{a numeric vector}
    \item{\code{difference}}{a numeric vector}
    \item{\code{year}}{a numeric vector}
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

\seealso{\code{\link{WillsDebates}}}

\examples{
data(WillsPresidents)
}

\keyword{datasets}
