\name{dataInfo}
\alias{dataInfo}
\docType{data}
\title{
\code{abd} Data Sets
}
\description{
Information about the location of data sets in 
\emph{Analysis of Biological Data}
}
\usage{data(dataInfo)}
\format{
  A data frame with 143 observations on the following 5 variables.
  \describe{
    \item{\code{name}}{name of data set}
    \item{\code{chapter}}{chapter in which data set appears}
    \item{\code{type}}{used in an \code{Example} or a \code{Problem}}
    \item{\code{number}}{example or problem number}
    \item{\code{sub}}{sub-problem:  \code{} \code{a} \code{b} \code{c}}
  }
}
\seealso{\code{\link{findData}}}

\examples{
data(dataInfo)
str(dataInfo)
}

\keyword{datasets}
