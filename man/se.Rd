\name{se}
\alias{se}
\title{Standard Error of the Mean}
\usage{
  se(x)
}
\arguments{
  \item{x}{a vector of numeric values}
}
\value{
  A single numeric value of the standard error of the mean
}
\description{
  \code{se} calculates the standard error of the mean of a
  vector.
}
\examples{
set.seed(2)
n <- 10
y <- rnorm(n)

sd(y)
sd(y)/sqrt(n)

se(y)
}
\author{
  Kevin M. Middleton (\email{kmm@csusb.edu})
}
\seealso{
  \code{\link{mean}} for the mean of a vector
  \code{\link{sd}} for the standard deviation of a vector
}
\keyword{univar}

