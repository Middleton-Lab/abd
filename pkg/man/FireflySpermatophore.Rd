\name{FireflySpermatophore}
\alias{FireflySpermatophore}
\docType{data}

\title{
Spermatophore Mass in Fireflies
}

\description{
Measurements of spermatophore mass (milligrams) in 35 fireflies 
(\emph{Photinus ignitus}).
}

\usage{data(FireflySpermatophore)}
\format{
A data frame with 35 observations of one variable.
\describe{
\item{\code{sp.mass}}{}
}
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Cratsley, C.K. and S.M. Lewis. 2003. Female preference for male courtship flashes in \emph{Photinus ignitus} fireflies. \emph{Behavioral Ecology} 14: 135-140.
}

\references{
\url{http://beheco.oxfordjournals.org/cgi/content/abstract/14/1/135}

\url{http://en.wikipedia.org/wiki/Firefly}
}

\examples{
data(FireflySpermatophore)
str(FireflySpermatophore)
histogram(~sp.mass, FireflySpermatophore, n=12)
}

\keyword{datasets}
