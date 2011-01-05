\name{sample}
\alias{sample}
\alias{sample.data.frame}
\alias{sample.default}
\docType{package}
\title{Random Samples and Permutations}
\description{
\code{sample} takes a sample of the specified size from the elements of 
\code{x} either with or without replacement.  Unlike \code{base::sample},
\code{sample} handles data frames in a reasonable manner.
}

\usage{
sample(x, size, replace = FALSE, prob = NULL)
}
\arguments{
  \item{x}{
  an integer, vector, or data frame.
}
  \item{size}{
a non-negative integer giving the number of items to choose.
}
  \item{replace}{
  logical indicating whether sampling should be performed with or without
  replacement.
}
  \item{prob}{
a vector of probability weights for obtaining the elements 
of the vector being sampled.
}
}



\value{
If \code{x} is a vector, a vector of length \code{size} 
with elements drawn from \code{x}.
If \code{x} is a integer, a vector of length \code{size} 
with elements drawn from \code{1:x}.
If \code{x} is a data frame, a data frame with \code{size} rows, all the columns
in \code{x} and one additional column indicating the location of the selected
rows in \code{x}.
}

\author{
Randall Pruim (\email{rpruim@calvin.edu})
}

\seealso{
\code{\link{sample}}
}
\examples{
x <- data.frame(letter=letters[1:10], number=1:10)
sample(x,3)
}
\keyword{manip}
