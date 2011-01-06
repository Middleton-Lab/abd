\name{Selection}
\alias{Selection}
\docType{data}

\title{
Data for Meta-analysis
}

\description{
Data for meta-analysis on the relationship between testosterone and aggression.
}

\usage{data(Selection)}
\format{
  A data frame with 814 observations on the following 8 variables.
  
  \describe{
	  \item{\code{species}}{species investigated}
	  \item{\code{traitname}}{trait investigated}
	  \item{\code{strength.of.selection}}{strength of selection}
	  \item{\code{sample.size}}{size of sample}
	  \item{\code{authors}}{authors of publication}
	  \item{\code{year}}{year of publication}
	  \item{\code{journal}}{journal of publication}
	  \item{\code{volume.pages}}{volume and pages}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Kingsolver, J.G., H.E. Hoekstra, J.M. Hoekstra, D. Berrigan, S.N. Vignieri, C.E. Hill, A. Hoang, P. Gibert, and P. Beerli. 2001. The strength of phenotypic selection in natural populations. \emph{The American Naturalist} 157: 245-261.}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(Selection)
histogram(~strength.of.selection, Selection,n=40)
table(Selection$species) -> s
table(s)
s[s>10] # most common species
table(Selection$traitname) -> t
table(t)
t[t>10] # most common traits
}

\keyword{datasets}
