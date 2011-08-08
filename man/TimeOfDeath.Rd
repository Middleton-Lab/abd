\name{TimeOfDeath}
\alias{TimeOfDeath}
\docType{data}

\title{
Hypoxanthine and Time Since Death
}

\description{
Hypoxanthine levels in the vitreous humour of the eye and time since death (hours) for 48 subjects.
}

\usage{data(TimeOfDeath)}
\format{
  A data frame with 48 observations on the following 2 variables.
  \describe{
    \item{\code{hours}}{a numeric vector}
    \item{\code{hypoxanthine}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
James, R.A., P.A. Hoadley, and B.G. Sampson. 1997. Determination of postmortem interval by sampling vitreous humor. \emph{American Journal of Forensic Medicine and Pathology} 18: 158-162.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(TimeOfDeath)
xyplot(hypoxanthine ~ hours, TimeOfDeath, type=c('p','r'))
}

\keyword{datasets}
