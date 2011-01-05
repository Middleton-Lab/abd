\name{WeddellSeals}
\alias{WeddellSeals}
\docType{data}

\title{
Energetic Cost of Diving
}

\description{
Comparison of oxygen consumption in feeding vs. non-feeding dives of the same length in the Weddell seal (\emph{Leptonychotes weddellii}).
}

\usage{data(WeddellSeals)}
\format{
  A data frame with 10 observations on the following 3 variables.
  \describe{
    \item{\code{individual}}{a numeric vector}
    \item{\code{oxygen.use.nonfeeding}}{a numeric vector}
    \item{\code{oxygen.use.feeding}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Williams, T.M., L.A. Fuiman, M. Horning, and R.W. Davis. 2004. The cost of foraging by a marine predator, the Weddell seal \emph{Leptonychotes weddellii}: pricing by the stroke. \emph{Journal of Experimental Biology} 207: 973-982.
}

\references{
\url{http://jeb.biologists.org/cgi/content/full/207/6/973}
}

\examples{
data(WeddellSeals)
xyplot(oxygen.use.nonfeeding ~ oxygen.use.feeding, WeddellSeals)
}

\keyword{datasets}
