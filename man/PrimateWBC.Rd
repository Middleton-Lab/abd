\name{PrimateWBC}
\alias{PrimateWBC}
\docType{data}

\title{
Primate White Blood Cell Counts and Promiscuity
}

\description{
White blood cell (WBC) counts in pairs of closely related primate species
}

\usage{data(PrimateWBC)}
\format{
  A data frame with 9 observations on the following 2 variables.
  \describe{
    \item{\code{WBC.less}}{a numeric vector}
    \item{\code{WBC.more}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Nunn, C.L., J.L. Gittleman, and J. Antonovics. 2000. Promiscuity and the primate immune system. \emph{Science} 290: 1168-1170.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(PrimateWBC)
xyplot(WBC.more ~ WBC.less, PrimateWBC)
}

\keyword{datasets}
