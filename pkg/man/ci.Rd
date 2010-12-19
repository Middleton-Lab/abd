\name{ci}

\alias{ci}
\alias{print.ci}

\title{
Confidence Interval for the Mean of a Normal Distribution
}
\description{
This function calculates the confidence interval (default is 95\%) for the mean of a normal distribution.
}
\usage{
  ci(x, conf.level = 0.95)
}

\arguments{
  \item{x}{a vector of values}
  \item{conf.level}{the confidence level for interval estimation; defaults to 0.95}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
  \item{x}{the original data \code{x}, for potential use}
  \item{x.bar}{the mean of \code{x}}
  \item{n}{the number of observations in \code{x}}
  \item{x.se}{the standard error of the mean of \code{x}}
  \item{df}{the degrees of freedom used to calculate the critical \emph{t}-statistic}
  \item{lower}{the lower bound of the confidence interval}
  \item{upper}{the upper bound of the confidence interval}
  \item{conf.level}{the confidence interval level}
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Kevin Middleton (\email{kmm@csusb.edu})
}

%\note{
%%  ~~further notes~~
%}

\seealso{
  \code{\link{mean}}, \code{\link{sd}}, \code{\link{se}}
}

\examples{
data(Stalkies)
ci(Stalkies)
}

\keyword{ univar }
