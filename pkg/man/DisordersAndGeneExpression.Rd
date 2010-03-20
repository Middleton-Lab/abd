\name{DisordersAndGeneExpression}
\alias{DisordersAndGeneExpression}
\docType{data}

\title{
Proteolipid Protein 1 Gene Expression
}

\description{
Expression levels of the proteolipid protein 1 gene (PLP1; \code{PLP1.expression}) in 45 individuals in one of three \code{group}s.
}

\usage{data(DisordersAndGeneExpression)}
\format{
  A data frame with 45 observations on the following 2 variables.
  \describe{
    \item{\code{group}}{a factor with levels: \code{bipolar}, \code{control}, and \code{schizo}}
    \item{\code{PLP1.expression}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Tkachev, D., M.L. Mimmack , M.M. Ryan , M. Wayland , T. Freeman , P.B. Jones , M. Starkey , M.J. Webster , R.H. Yolken , S. Bahn. 2003. Oligodendrocyte dysfunction in schizophrenia and bipolar disorder. \emph{Lancet} 362(9386): 798-805.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(DisordersAndGeneExpression)
str(DisordersAndGeneExpression)
}

\keyword{datasets}
