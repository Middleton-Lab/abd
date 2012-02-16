\name{as.xtabs}
\alias{as.xtabs}
\alias{as.xtabs.data.frame}
\alias{as.xtabs.matrix}
\title{Convert objects to xtabs format}
\usage{
  as.xtabs(x, ...)

  \method{as.xtabs}{data.frame} (x, rowvar = NULL,
    colvar = NULL, labels = 1, ...)

  \method{as.xtabs}{matrix} (x, rowvar = NULL,
    colvar = NULL, ...)
}
\arguments{
  \item{x}{object (typically a data frame) to be converted
  to \code{xtabs} format}

  \item{...}{additional arguments to be passed to or from
  methods.}

  \item{rowvar}{name of the row variable as character
  string}

  \item{colvar}{name of the column variable as character
  string}

  \item{labels}{Column of data frame that contains the
  lables of the row variable.}
}
\value{
  An \code{xtabs} object.
}
\description{
  Currently able to convert a data frame or a matrix into
  an \code{xtabs} object.
}
\details{
  The intended use it to convert a two-way contingency
  table stored in a data frame or a matrix into an
  \code{xtabs} object.
}
\examples{
# example from example(fisher.test)
df <- data.frame( X=c('Tea','Milk'), Tea=c(3,1), Milk=c(1,3) )
xt <- as.xtabs(df, rowvar="Guess", colvar="Truth"); xt
if (require(vcd)) { mosaic(xt) }
}
\author{
  Randall Pruim (\email{rpruim@calvin.edu})
}
\keyword{manip}

