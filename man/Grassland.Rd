\name{Grassland}
\alias{Grassland}
\docType{data}

\title{
Grassland Diversity
}

\description{
Species diversity in 10 experimental plots in the Park Grass Experiment at Rothamsted Experimental Station to which varying numbers of nutrients have been added.
}

\usage{data(Grassland)}
\format{
  A data frame with 10 observations on the following 2 variables.
  \describe{
    \item{\code{nutrients}}{a numeric vector}
    \item{\code{num.species}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{ Harpole, W. S. and D. Tilman. 2007. Grassland species loss due to
reduced niche dimension. \emph{Nature} 446: 791-793.  }

\references{
\url{http://www.rothamsted.ac.uk/}
}

\examples{
data(Grassland)
xyplot(num.species ~ jitter(nutrients, amount=0.1), Grassland, pch=16)
}

\keyword{datasets}
