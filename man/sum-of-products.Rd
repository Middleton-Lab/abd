\name{sum_of_squares}
\alias{sum_of_products}
\alias{sum_of_squares}
\title{Sum of Squares and Sum of Products}
\usage{
  sum_of_squares(x)

  sum_of_products(x, y)
}
\arguments{
  \item{x}{a vector of numeric values}

  \item{y}{a vector of numeric values}
}
\value{
  a numeric value
}
\description{
  Calculate the sum of squares (\code{sum_of_squares}) or
  sum of products (\code{sum_of_products}) of a numeric
  vector or pair of numeric vectors
}
\examples{
set.seed(4)
x <- rnorm(10)
sum_of_squares(x)

y <- rnorm(10)
sum_of_products(x, y)
}
\author{
  Kevin M. Middleton (\email{kmm@csusb.edu})
}
\keyword{univar}

