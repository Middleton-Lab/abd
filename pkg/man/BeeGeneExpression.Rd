\name{BeeGeneExpression}
\alias{BeeGeneExpression}
\docType{data}

\title{
Foraging Gene Expression
}

\description{
Levels of expression of the foraging gene (\emph{for}; \code{Expression}) in two worker types (\code{type}) in three bee colonies (\code{colony}). Note that \code{colony} is not coded as a factor.
}

\usage{data(BeeGeneExpression)}
\format{
  A data frame with 6 observations on the following 3 variables.
  \describe{
    \item{\code{type}}{a factor with levels \code{forager} \code{nurse}}
    \item{\code{colony}}{a numeric identifier}
    \item{\code{expression}}{expression level of the \emph{for} gene}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Ben-Shahar, Y., A. Robichon, M.B. Sokolowski, and G.E. Robinson. 2002. Influence of gene action across different time scales on behavior. \emph{Science} 296: 741-744.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(BeeGeneExpression)
str(BeeGeneExpression)
BeeGeneExpression
}

\keyword{datasets}
