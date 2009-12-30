\name{NorthSeaCodRecruits}
\alias{NorthSeaCodRecruits}
\docType{data}

\title{
Atlantic cod recruits
}

\description{
Size (log10 transformed; \code{log10.recruits}) of Atlantic cod (\emph{Gadus morhua}) recruits in the North Sea in 39 years
}

\usage{data(NorthSeaCodRecruits)}
\format{
  A data frame with 39 observations on the following variable.
  \describe{
    \item{\code{log10.recruits}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Beaugrand, G., K.M. Brander, J.A. Lindley, S. Souissi, and P.C. Reid. 2003. Plankton effect on cod recruitment in the North Sea. \emph{Nature} 426: 661-664.
}

\references{
\url{http://www.nature.com/nature/journal/v426/n6967/abs/nature02164.html}
}

\examples{
data(NorthSeaCodRecruits)
NorthSeaCodRecruits$log10.recruits
}

\keyword{datasets}
