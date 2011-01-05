\name{NorthSeaCod}
\alias{NorthSeaCod}
\docType{data}

\title{
Atlantic Cod Recruits
}

\description{
Number (\eqn{\log_{10}}{log10} transformed) of Atlantic cod (\emph{Gadus
morhua}) that recruited (grew to catchable size) in the North Sea over a 39
years span.
}

\usage{data(NorthSeaCod)}
\format{
A data frame with 39 observations of one variable.
\describe{
\item{\code{log10.recruits}}{}
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
data(NorthSeaCod)
favstats(NorthSeaCod$log10.recruits)
}

\keyword{datasets}
