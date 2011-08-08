\name{LionNoses}
\alias{LionNoses}
\docType{data}

\title{
Lion Age and Nose Coloration
}

\description{ Ages (in years; \code{age}) of 32 male lions and relative
coloration of their noses (\code{proportion.black}).  }

\usage{data(LionNoses)}
\format{
  A data frame with 32 observations on the following 2 variables.
  \describe{
    \item{\code{age}}{a numeric vector}
    \item{\code{proportion.black}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Whitman, K., A.M. Starfield, H.S. Quadling and C. Packer. 2004. Sustainable trophy hunting of African lions. \emph{Nature} 428: 175-178.
}

\references{
\url{http://www.nature.com/nature/journal/v428/n6979/abs/nature02395.html}
}

\examples{
data(LionNoses)
xyplot(age ~ proportion.black, LionNoses)
}

\keyword{datasets}
