\name{col.abd}
\alias{col.abd}
\title{
Lattice theme for Analysis of Biological Data
}
\description{
This theme will help produce plots with color scheme similar to the one used in
\emph{Analysis of Biological Data}
}
\usage{
col.abd(bw = FALSE, lty = 1:7)
}
\arguments{
  \item{bw}{a logical. Use a grayscale theme instead of color?}
  \item{lty}{line types used for \code{panel.superpose}}
}
\value{
a list that can be used as a \code{lattice} theme.
}
\author{
Randall Pruim (\email{rpruim@calvin.edu})
}

\seealso{
\code{\link{trellis.par.set}}
\code{\link{show.settings}}
\code{\link{col.whitebg}}
}
\examples{
trellis.par.set(theme=col.abd(bw=TRUE))
show.settings()
trellis.par.set(theme=col.abd(lty=1))
show.settings()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphics }
