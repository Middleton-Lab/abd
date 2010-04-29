\name{CocaineHighVolkowEtal}
\alias{CocaineHighVolkowEtal}
\docType{data}

\title{
Effects of Cocaine on Dopamine Receptors
}

\description{
Percent of dopamine receptors blocked (\code{percent.bocked}) and the perceived level of high as determined by PET scans (\code{high}) in 34 humans.
}

\usage{data(CocaineHighVolkowEtal)}
\format{
  A data frame with 34 observations on the following 2 variables.
  \describe{
    \item{\code{percent.blocked}}{a numeric vector}
    \item{\code{high}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Volkow, N.D., G.-J. Wang, R.W. Foltin, J.S. Fowler, N.N. Abumrad,  S. Vitkun, J. Logan, S.J. Gatley, N. Pappas, R. Hitzemann, and C.E. Shea. 1997. Relationship between subjective effects of cocaine and dopamine transporter occupancy. \emph{Nature} 386: 827-830.
}

\references{
\url{http://www.nature.com/nature/journal/v386/n6627/abs/386827a0.html}
}

\examples{
data(CocaineHighVolkowEtal)
str(CocaineHighVolkowEtal)
CocaineHighVolkowEtal
}

\keyword{datasets}
