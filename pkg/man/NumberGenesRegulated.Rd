\name{NumberGenesRegulated}
\alias{NumberGenesRegulated}
\docType{data}

\title{
NumberGenesRegulated
}

\description{
Number of genes regulated by 109 regulatory genes of \emph{Saccharomyces
  cerevisiae}.
}

\usage{data(NumberGenesRegulated)}
\format{
  A data frame with 26 observations on the following 2 variables.
  \describe{
    \item{\code{ngenes}}{a numeric vector}
    \item{\code{frequency}}{a numeric vector}
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
data(NumberGenesRegulated)
str(NumberGenesRegulated)
NumberGenesRegulated
}

\keyword{datasets}
