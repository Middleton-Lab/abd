\name{se}
\alias{se}
\title{
Standard Error of the Mean
}

\description{
Calculate the standard error of the mean of a vector.
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
A single value of the standard error of the mean
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Kevin Middleton
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mean}}, \code{\link{sd}}
}

\examples{
set.seed(2)
n <- 10
y <- rnorm(n)
mean(y)
sd(y)
sd(y)/sqrt(n)
se(y)
}

\keyword{ univar }
