\name{histochart}
\alias{histochart}
\title{
	Histogram from tabulated data
}
\description{
Uses \code{\link{barchart}} to build a histogram from tabulated data.
}
\usage{
histochart(x, data=NULL, box.ratio = 100, origin = 0, horizontal = FALSE, ...)
}
\arguments{
  \item{x}{ formula of form \code{frequency ~ value}
}
  \item{data}{ data frame in which the formula \code{x} is interpreted
}
  \item{box.ratio}{ ratio of bar widths to gaps between bars
}
  \item{origin}{
  where do bars begin?
}
  \item{horizontal}{ Should bars go horizontal?
}
  \item{\dots}{
  other arguments passed to \code{\link{barchart}}
}
}
\details{
This is just a convenience wrapper around \code{\link{barchart}}.
}
\author{
Randall Pruim (\email{rpruim@calvin.edu})
}

\seealso{
\code{\link{barchart}}
}
\examples{
histochart( dbinom(0:30, 30, 0.35) ~ 0:30 )
}
\keyword{ graphics }
