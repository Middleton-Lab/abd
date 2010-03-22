\name{sum_of_squares}

\alias{sum_of_squares}
\alias{sum_of_products}

\title{
Sum of Squares and Sum of Products
}

\description{
Calculate the sum of squares (\code{sum_of_squares}) or sum of products (\code{sum_of_products}) of a numeric vector or pair of numeric vectors
}

\usage{
sum_of_squares(x)

sum_of_products(x, y)
}

\arguments{
  \item{x}{a vector of numeric values}
  \item{y}{a vector of numeric values}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
a numeric value
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

%\seealso{
%\code{\link{sum_of_products}}
%}

\examples{
set.seed(4)
x <- rnorm(10)
sum_of_squares(x)

y <- rnorm(10)
sum_of_products(x, y)
}

\keyword{ univar }
