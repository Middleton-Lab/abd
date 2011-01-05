\name{Lobsters}
\alias{Lobsters}
\docType{data}

\title{
Lobster Orientation
}

\description{
Orientation of 15 lobsters relative to initial position.
}

\usage{data(Lobsters)}
\format{
A data frame with 15 observations of one variable.
\describe{
\item{\code{orientation}}{}
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Boles, L.C. and K.J. Lohmann. 2003. True navigation and magnetic maps in spiny lobsters. \emph{Nature} 421: 60-63.
}

\references{
\url{http://www.unc.edu/depts/geomag/PDFGeomag/BolesandLohmann2003.pdf}
}

\examples{
data(Lobsters)
histogram(~orientation, Lobsters)
dotplot(~orientation, Lobsters)
}

\keyword{datasets}
