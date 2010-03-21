\name{sum.of.squares}

\alias{sum.of.squares}
\alias{sum.of.products}

\title{
Sum of Squares and Sum of Products
}

\description{
Calculate the sum of squares (\code{sum.of.squares}) or sum of products (\code{sum.of.products}) of a numeric vector or pair of numeric vectors
}

\usage{
sum.of.squares(x)

sum.of.products(x, y)
}

\arguments{
  \item{x}{a vector of values}
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

\seealso{
\code{\link{sum.of.products}}
}

\examples{
set.seed(4)
x <- rnorm(10)
sum.of.squares(x)

y <- rnorm(10)
sum.of.products(x, y)
}

\keyword{ univar }
