\name{cumfreq}
\alias{cumfreq}
\alias{cumfreq.default}
\alias{cumfreq.formula}
\alias{panel.cumfreq}
\alias{prepanel.cumfreq}
\title{Cumulative frequency plots}
\usage{
  cumfreq(x, data, ...)

  \method{cumfreq}{formula} (x, data = NULL, subscripts,
    ...)

  \method{cumfreq}{default} (x, ...)

  prepanel.cumfreq(x, ...)

  panel.cumfreq(x, type = c("smooth", "step"),
    groups = NULL, ...)
}
\arguments{
  \item{x}{a formula or numeric vector}

  \item{data}{a data frame in which \code{x} is evaluated
  if \code{x} is a formula.}

  \item{\dots}{other lattice arguments}

  \item{subscripts}{as in lattice plots}

  \item{type}{smooth or step-function?}

  \item{groups}{grouping variable}
}
\description{
  A high-level function for producing a cumulative
  frequency plot using \code{lattice} graphics.
}
\examples{
cumfreq(~count, DesertBirds, xlab = 'Species Abundance')
}
\author{
  Randall Pruim (\email{rpruim@calvin.edu})
}
\seealso{
  \code{\link{histogram}}, \code{\link{densityplot}}
}
\keyword{graphics}

