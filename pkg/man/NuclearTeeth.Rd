\name{NuclearTeeth}
\alias{NuclearTeeth}
\docType{data}

\title{
Radioactive Teeth
}

\description{
Actual birth year and birth year estimated from relative radioactivity of the enamel for 20 samples.
}

\usage{data(NuclearTeeth)}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{actual}}{a numeric vector}
    \item{\code{estimated}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Spalding, K.L., B.A. Buchholz, L.-E. Bergman, H. Druid, and J. Frisén. 2005. Age written in teeth by nuclear tests. \emph{Nature} 437: 333-334.
}

\references{
\url{http://www.nature.com/nature/journal/v437/n7057/full/437333a.html}
}

\examples{
data(NuclearTeeth)
str(NuclearTeeth)
NuclearTeeth
}

\keyword{datasets}
