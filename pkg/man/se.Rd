\name{se}
\alias{se}
\title{
Standard Error of the Mean
}

\description{
\code{se} calculates the standard error of the mean of a vector.
}

\usage{
se(x)
}

\arguments{
  \item{x}{a vector of numeric values}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A single numeric value of the standard error of the mean
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Kevin M. Middleton (kmm@csusb.edu)
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mean}} for the mean of a vector
\code{\link{sd}} for the standard deviation of a vector
}

\examples{
set.seed(2)
n <- 10
y <- rnorm(n)

sd(y)
sd(y)/sqrt(n)

se(y)
}

\keyword{ univar }
