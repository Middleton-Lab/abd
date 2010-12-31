\name{cumfreq}
\alias{cumfreq}
\alias{cumfreq.formula}
\alias{cumfreq.default}
\alias{panel.cumfreq}
\alias{prepanel.cumfreq}
\title{
Cumulative Frequency Plots
}
\description{
High-level graphics functions for creating cumulative frequency plots
using \code{lattice} graphics.
}
\usage{
cumfreq(x, data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector or a formula }
  \item{data}{ if \code{x} is a formula, then \code{data} is a 
  data frame in which the formula is interpreted.
}
  \item{\dots}{
  other arguments passed to \code{panel.cumfreq}.
}
}
\value{
A trellis object representing the cumulative frequency plot.
}
\author{ Randall Pruim }

\examples{
cumfreq(~Sepal.Length, groups=Species, data=iris)
cumfreq(~Sepal.Length, groups=Species, data=iris, type='step')
}
\keyword{ graphics }
