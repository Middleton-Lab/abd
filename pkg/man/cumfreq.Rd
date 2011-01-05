\name{cumfreq}
\alias{cumfreq}
\alias{cumfreq.default}
\alias{panel.cumfreq}
\alias{prepanel.cumfreq}
\alias{cumfreq.formula}
\title{
Cumulative frequency plots
}
\description{
A high-level function for producing a cumulative frequency plot using 
\code{lattice} graphics.
}
\usage{
cumfreq(x, data, ...)
}
\arguments{
  \item{x}{ a formula or numeric vector }
  \item{data}{ a data frame in which \code{x} is evaluated if \code{x} is a formula.
}
  \item{\dots}{
  other lattice arguments
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\author{
Randall Pruim (\email{rpruim@calvin.edu})
}


\seealso{
\code{\link{histogram}}, \code{\link{densityplot}}
}
\examples{
cumfreq(~count,DesertBirds, xlab='Species Abundance')
}
\keyword{ graphics }
