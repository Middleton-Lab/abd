\name{findData}
\alias{findData}
\title{
Find data in \emph{Analysis of Biological Data}
}
\description{
A utility function to assist users to locate data sets in
\emph{Analysis of Biological Data} within the \code{abd} package.
}
\usage{
findData(..., 
	chapters = 1:21, 
	types = c("Example", "Problem"), 
	numbers = 1:100, 
	pattern = "*", 
	ignore.case = TRUE)
}
\arguments{
  \item{\dots}{
  values for any of 
  \code{chapters},
  \code{types}, or 
  \code{pattern}.
  Which is meant will be inferred from the type of object supplied.
  This allows users to specify these values in any order and without
  naming.
}
  \item{chapters}{
  		a numeric vector of chapters to search within
}
  \item{types}{
  		a sub-vector of c('Example','Problem')
}
  \item{numbers}{
  		a numeric vector of problem numbers
}
  \item{pattern}{
  		a pattern to use for regular expression matching against the name
		of the data frame.
}
  \item{ignore.case}{
  		should case be ignored when matching pattern?
}
}
\value{
A data frame describing data sets from \code{abd} that match the search criteria,
or NULL if there are no matches.
}
\author{
Randall Pruim (\email{rpruim@calvin.edu})
}

\seealso{
\code{\link{dataInfo}}
}
\examples{
# find all data from examples in chapters 3 and 4
findData(3:4, 'Example')
# order doesn't matter
findData('Example', 3:4)
# look for data sets with Example in their name.
findData(pattern='Example')
# look for data sets with Exercise in their name.
findData('Exercise')
}
\keyword{ datasets }
