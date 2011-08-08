\name{MoleRats}
\alias{MoleRats}
\docType{data}

\title{
Energy Expenditure in Mole Rats
}

\description{
Energy expenditure (\code{ln.energy}) in two castes (\code{caste}) of Damaraland mole rats (\emph{Cryptomys damarensis}) with body mass (\code{ln.mass}) as a covariate.
}

\usage{data(MoleRats)}
\format{
  A data frame with 35 observations on the following 3 variables.
  \describe{
    \item{\code{caste}}{a factor with levels \code{lazy} and \code{worker}}
    \item{\code{ln.mass}}{a numeric vector}
    \item{\code{ln.energy}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Scantlebury, M., J.R. Speakman, M.K. Oosthuizen, T.J. Roper and N.C. Bennett. 2006. Energetics reveals physiologically distinct castes in a eusocial mammal. \emph{Nature} 440: 795-797.
}

\references{
\url{http://www.nature.com/nature/journal/v440/n7085/abs/nature04578.html}
}

\examples{
data(MoleRats)
MoleRats
}

\keyword{datasets}
