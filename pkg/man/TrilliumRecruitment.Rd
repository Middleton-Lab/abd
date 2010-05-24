\name{TrilliumRecruitment}
\alias{TrilliumRecruitment}
\docType{data}

\title{
Trillium Recruitment near Clearcuts
}

\description{
Recruitment of \emph{Trillium} and distance to nearest clearcut in eight populations in southwestern Oregon.
}

\usage{data(TrilliumRecruitment)}
\format{
  A data frame with 8 observations on the following 3 variables.
  \describe{
    \item{\code{population}}{a numeric vector}
    \item{\code{edgedist}}{a numeric vector}
    \item{\code{recruitment}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Jules, E.S. and B.J. Rathcke. 1999.  Mechanisms of reduced trillium recruitment along edges of old-growth forest fragments. \emph{Conservation Biology} 13:784-793
}

\references{
\url{http://www.jstor.org/stable/2641693}
}

\examples{
data(TrilliumRecruitment)
str(TrilliumRecruitment)
TrilliumRecruitment
}

\keyword{datasets}
