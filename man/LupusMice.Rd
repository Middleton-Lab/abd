\name{LupusMice}
\alias{LupusMice}
\docType{data}

\title{
Autoimmune Reactivity in Lupus-prone Mice
}

\description{
Autoimmune reactivity (\code{dilution} at which reactivity could be detected) in three \code{treatment}s of lupus-prone mice.
}

\usage{data(LupusMice)}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels: \code{enhanced}, \code{sham}, and \code{untreated}}
    \item{\code{dilution}}{a numeric vector of the dilution level at which reactivity could be detected}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
McGaha, T.L., B. Sorrentino, and J.V. Ravetch. 2005. Restoration of tolerance in lupus by targeted inhibitory receptor expression. \emph{Science} 307: 590-593.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(LupusMice)
str(LupusMice)
}

\keyword{datasets}
