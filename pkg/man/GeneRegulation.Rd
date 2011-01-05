\name{GeneRegulation}
\alias{GeneRegulation}
\docType{data}

\title{
Gene Regulation in Saccharomyces
}

\description{
Number of genes regulated by 109 regulatory genes of \emph{Saccharomyces
  cerevisiae}.
}

\usage{data(GeneRegulation)}
\format{
  A data frame with 26 observations on the following 2 variables.
  \describe{
    \item{\code{genes.regulated}}{a numeric vector}
    \item{\code{count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Guelzim, N., S. Bottani, P. Bourgine and F. Képès. 2002. Topological and
  causal structure of the yeast transcriptional regulatory network.
  \emph{Nature Genetics} 31: 60-63.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(GeneRegulation)
str(GeneRegulation)
xyplot(count ~ genes.regulated, GeneRegulation, type='h', lwd=3)
}

\keyword{datasets}
