\docType{data}
\name{Greenspace}
\alias{Greenspace}
\title{Diversity in Urban Green Space}
\format{A data frame with 15 observations on the following 6 variables.
\describe{ \item{site}{a factor with levels \code{A} - \code{O}}
\item{attachment}{a numeric vector} \item{area}{a numeric
vector} \item{butterfly}{a numeric vector} \item{bird}{a
numeric vector} \item{ln.plant}{a numeric vector} }}
\source{
  Fuller, R.A., K.N. Irvine, P. Devine-Wright, P.H. Warren,
  and K.J. Gaston. 2007. Psychological benefits of
  greenspace increase with biodiversity. \emph{Biology
  Letters} 3: 390-394.
}
\description{
  Measures of biodiversity in 15 urban green spaces in
  Sheffield, England.
}
\examples{
data(Greenspace)
str(Greenspace)
splom(Greenspace[,2:6])
}
\references{
  \url{http://rsbl.royalsocietypublishing.org/content/3/4/390.abstract}
}
\keyword{datasets}

