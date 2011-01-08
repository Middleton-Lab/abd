\name{Hemoglobin}
\alias{Hemoglobin}
\docType{data}

\title{
Hemoglobin Levels in High Altitude Populations
}

\description{
Relative rates of hemoglobin concentration in four populations of humans living at different altitudes.
}

\usage{data(Hemoglobin)}
\format{
  A data frame with 40 observations on the following 3 variables.
  \describe{
    \item{\code{hemoglobin}}{a numeric vector}
    \item{\code{group}}{a factor with levels: \code{Andes}, \code{Ethiopia}, \code{Tibet}, and \code{USA}}
    \item{\code{relative.frequency}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Beall, C.M., M.J. Decker, G.M. Bittenham, I. Kushner, A. Gebremedhin, K.P. Strohl. 2002. An Ethiopian pattern of human adaptation to high-altitude hypoxia. \emph{Proceeding of the National Academy of Sciences (USA)} 99(26): 17215-17218.}

\references{
\url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC139295/}
}

\examples{
data(Hemoglobin)
str(Hemoglobin)

xyplot(relative.frequency ~ hemoglobin | group, Hemoglobin,
	type ='h', lwd=4, layout=c(1,4))
}

\keyword{datasets}
