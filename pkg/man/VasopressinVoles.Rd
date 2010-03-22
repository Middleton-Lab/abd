\name{VasopressinVoles}
\alias{VasopressinVoles}
\docType{data}

\title{
Vasopressin Manipulation in the Meadow Vole
}

\description{
Time spent with a female (\code{percent}) in control and vasopressin-enhanced groups (\code{treatment}) of meadow voles (\emph{Microtus pennsylvanicus}).
}

\usage{data(VasopressinVoles)}
\format{
  A data frame with 31 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{control} and \code{enhanced}}
    \item{\code{percent}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
\emph{inferred from} Lim, M.M., Z. Wang, D.E. \enc{Olazábal}{Olazabal}, X. Ren, E.F. Terwilliger, and L.J. Young. 2004. Enhanced partner preference in a promiscuous species by manipulating the expression of a single gene. \emph{Nature} 429: 754-757.
}

\references{
\url{http://www.nature.com/nature/journal/v429/n6993/full/nature02539.html}
}

\examples{
data(VasopressinVoles)
VasopressinVoles
}

\keyword{datasets}
