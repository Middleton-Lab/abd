\name{RopeTrick}
\alias{RopeTrick}
\docType{data}

\title{
Indian Rope Trick
}

\description{
Perceived impressiveness (\code{impressiveness}) of a written account of
the Indian Rope Trick and the corresponding number of \code{year}s since it was
witnessed.
}

\usage{data(RopeTrick)}
\format{
  A data frame with 21 observations on the following 2 variables.
  \describe{
    \item{\code{years}}{a numeric vector}
    \item{\code{impressiveness}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Wiseman, R. and P. Lamont. 1996. Unravelling the Indian rope-trick. \emph{Nature} 383: 212-213.
}

\references{
\url{http://www.richardwiseman.com/resources/ropeJSPR.pdf}
}

\examples{
data(RopeTrick)
xyplot(impressiveness ~ years, RopeTrick)

}

\keyword{datasets}
