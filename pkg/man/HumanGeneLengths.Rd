\name{HumanGeneLengths}
\alias{HumanGeneLengths}
\docType{data}

\title{
Human Gene Lengths
}

\description{
Lengths in number of nucleotides (\code{gene.length}) for 20,290 human genes
}

\usage{data(HumanGeneLengths)}
\format{
  A data frame with 20,290 observations on the following variable.
  \describe{
    \item{\code{gene.length}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Hubbard, T., D. Andrews, M. Caccamo, G. Cameron, Y. Chen, M. Clamp, L. Clarke, G. Coates, T. Cox, F. Cunningham, V. Curwen, T. Cutts, T. Down, R. Durbin, X. M. Fernandez-Suarez, J. Gilbert, M. Hammond, J. Herrero, H. Hotz, K. Howe, V. Iyer, K. Jekosch, A. Kahari, A. Kasprzyk, D. Keefe, S. Keenan, F. Kokocinsci, D. London, I. Longden, G. McVicker, C. Melsopp, P. Meidl, S. Potter, G. Proctor, M. Rae, D. Rios, M. Schuster, S. Searle, J. Severin, G. Slater, D. Smedley, J. Smith, W. Spooner, A. Stabenau, J. Stalker, R. Storey, S. Trevanion, A. Ureta-Vidal, J. Vogel, S. White, C. Woodwark, and E. Birne. 2005. Ensembl 2005. \emph{Nucleic Acids Research} 33: D447-D453.
}

\references{
\url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC540092/}

\url{http://www.ensembl.org/}
}

\examples{
data(HumanGeneLengths)
str(HumanGeneLengths)
histogram(~gene.length, HumanGeneLengths, subset=gene.length<15000)
}

\keyword{datasets}
