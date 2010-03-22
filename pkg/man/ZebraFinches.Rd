\name{ZebraFinches}
\alias{ZebraFinches}
\docType{data}

\title{
Zebra Finch Carotenoids
}

\description{
Data on cell-mediated immunocompetence (\code{PHA}) and humoral immunity (\code{SRBC}) in Zebra Finches that received supplemental carotenoids (\code{CAROT}) and those that did not (\code{NO}).
}

\usage{data(ZebraFinches)}
\format{
  A data frame with 20 observations on the following 3 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels: \code{CAROT} and \code{NO}}
    \item{\code{PHA}}{a numeric vector}
    \item{\code{SRBC}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
McGraw, K.J. and D.R. Ardia. 2003. Carotenoids, immunocompetence, and the information content of sexual colors: an experimental test. \emph{The American Naturalist} 162: 704-712.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(ZebraFinches)
ZebraFinches
}

\keyword{datasets}
