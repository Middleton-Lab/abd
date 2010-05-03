\name{odds.ratio}
\alias{odds.ratio}
\alias{print.odds.ratio}

\title{
Odds Ratio for 2X2 Contingency Tables
}
\description{
This function calculates the odds ratio for a 2 X 2 contingency table and a confidence interval for the estimated odds ratio. \code{x} should be a matrix or data.frame. "Successes" should be located in column 1 of \code{x}, and the treatment of interest should be located in row 2. The odds ratio is calculated as (Odds row 2) / (Odds row 1). The confidence interval is calculated from the log(OR) and backtransformed.
}

\usage{
odds.ratio(x, conf.int = 0.95)
}

\arguments{
  \item{x}{a 2 X 2 matrix or data.frame of counts}
  \item{conf.int}{the confidence interval level}
}

%\details{
%
%}

\value{
  \item{p1, p2}{Proportions for rows 1 and 2}
  \item{o1, o2}{Odds for rows 1 and 2}
  \item{OR}{Odds ratio}
  \item{lower}{the lower bound of the confidence interval}
  \item{upper}{the upper bound of the confidence interval}
  \item{conf.int}{the confidence interval level}
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Kevin Middleton (kmm@csusb.edu)
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{chisq.test}}
}

\examples{
M1 <- matrix(c(14, 38, 51, 11), nrow = 2)
M1
odds.ratio(M1)

M2 <- matrix(c(18515, 18496, 1427, 1438), nrow = 2)
M2
odds.ratio(M2)
}

\keyword{ univar }
