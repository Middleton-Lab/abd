\name{SockeyeFemales}
\alias{SockeyeFemales}
\docType{data}

\title{
Body Masses of Female Sockeye Salmon
}

\description{ Body Masses of 228 female Sockeye Salmon 
(\emph{Oncorhynchus nerka}; 
\url{http://www.nmfs.noaa.gov/pr/species/fish/sockeyesalmon.htm}) 
}

\usage{data(SockeyeFemales)}
\format{A data frame with 228 observations of a single variable.
\describe{
\item{\code{mass}}{body mass (kg)}
}
}


\source{
Hendry, A.P., O.K. Berg, and T.P. Quinn. 1999. Condition dependence and adaptation-by-time: Breeding date, life history, and energy allocation within a population of salmon. \emph{Oikos} 85: 499-514.
}

\references{
\url{http://www.jstor.org/stable/3546699}
}

\examples{
data(SockeyeFemales)
str(SockeyeFemales)
summary(SockeyeFemales)
}

\keyword{datasets}
