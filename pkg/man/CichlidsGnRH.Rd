\name{CichlidsGnRH}
\alias{CichlidsGnRH}
\docType{data}

\title{
GnRH Levels in Cichlids
}

\description{ Levels of mRNA for gonadotropin-releasing hormone in cichlids
(\emph{Haplochromis burtoni}) that are (\emph{n} = 5) and are not
(\emph{n} = 6) territorial.  }

\usage{data(CichlidsGnRH)}
\format{
  A data frame with 11 observations on the following 2 variables.
  \describe{
    \item{\code{territorial}}{a factor with levels \code{No} and \code{Yes}}
    \item{\code{GnRH.mRNA}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
White, S.A., T. Nguyen, and R.D. Fernald. 2002. Social regulation of gonadotropin-releasing hormone. \emph{Journal of Experimental Biology} 205: 2567-2581.
}

\references{
\url{http://jeb.biologists.org/cgi/content/abstract/205/17/2567}
}

\examples{
data(CichlidsGnRH)
xyplot(GnRH.mRNA ~ territorial, CichlidsGnRH, type=c('p','a'))
}

\keyword{datasets}
