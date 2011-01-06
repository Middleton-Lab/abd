\name{YeastGenes}
\alias{YeastGenes}
\docType{data}

\title{
Yeast Regulatory Genes
}

\description{
Number of genes regulated by 109 yeast regulatory genes.
}

\usage{data(YeastGenes)}
\format{
  A data frame with 6 observations on the following 2 variables.
  \describe{
    \item{\code{genes.controlled}}{a numeric vector}
    \item{\code{count}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Guelzim, N., S. Bottani, P. Bourgine and F. Képès. 2002. Topological and causal structure of the yeast transcriptional regulatory network. \emph{Nature Genetics} 31: 60-63.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(YeastGenes)
str(YeastGenes)
barchart(count ~ genes.controlled , YeastGenes, horizontal=FALSE)
}

\keyword{datasets}
