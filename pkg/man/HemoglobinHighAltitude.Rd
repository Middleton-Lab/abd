\name{HemoglobinHighAltitude}
\alias{HemoglobinHighAltitude}
\docType{data}

\title{
HemoglobinHighAltitude
}

\description{
Relative rates of hemoglobin concentration in four populations of humans.
}

\usage{data(HemoglobinHighAltitude)}
\format{
  A data frame with 40 observations on the following 3 variables.
  \describe{
    \item{\code{haemoglobin}}{a numeric vector}
    \item{\code{group}}{a factor with levels \code{Andes}, \code{Ethiopia}, \code{Tibet}, and \code{USA}}
    \item{\code{relative.frequency}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Beall, C.M., M.J. Decker, G.M. Bittenham, I. Kushner, A. Gebremedhin, K.P. Strohl. 2002. An Ethiopian pattern of human adaptation to high-altitude hypoxia. \emph{Proc Natl Acad Sci U S A} 99(26): 17215-17218.}

\references{
\url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC139295/}
}

\examples{
data(HemoglobinHighAltitude)

str(HemoglobinHighAltitude)

# Fig. 2.4-1
require(ggplot2)
p <- ggplot(HemoglobinHighAltitude, aes(haemoglobin, relative.frequency))
p + geom_bar(stat="identity", fill = 'red')  +
  facet_grid(group ~ .) +
  scale_x_continuous('Hemoglobin concentration (g/dL)') +
  scale_y_continuous('Relative frequency')

# TODO
# Fig. 2.4-2

}

\keyword{datasets}
