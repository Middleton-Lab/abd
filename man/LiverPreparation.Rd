% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{LiverPreparation}
\alias{LiverPreparation}
\title{Liver Preparation}
\format{
A data frame with 5 observations on the following 2 variables.
\describe{ \item{concentration}{a numeric vector}
\item{unbound.fraction}{a numeric vector} }
}
\source{
Smallwood, R.H., D.J. Morgan, G.W. Mihaly, and R.A. Smallwood. 1998.
Effect of a protein binding change on unbound and total plasma
concentrations for drugs of intermediate hepatic extraction. \emph{Journal
of Pharmacokinetics and Pharmacodynamics} 16:397-411.
}
\usage{
LiverPreparation
}
\description{
The unbound fraction of taurocholate for each of five concentrations of
administered taurocholate.
}
\examples{

str(LiverPreparation)
xyplot(unbound.fraction ~ concentration, LiverPreparation)


}
\references{
\url{http://www.ncbi.nlm.nih.gov/pubmed/3199317}
}
\keyword{datasets}
