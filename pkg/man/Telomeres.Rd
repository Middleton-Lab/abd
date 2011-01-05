\name{Telomeres}
\alias{Telomeres}
\docType{data}

\title{
Telomere Shortening
}

\description{
Telomere length (ratio) and years since their child's diagnosis with chronic illness.
}

\usage{data(Telomeres)}
\format{
  A data frame with 39 observations on the following 2 variables.
  \describe{
    \item{\code{years}}{a numeric vector}
    \item{\code{telomere.length}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Epel, E.S., E.H. Blackburn, J. Lin, F.S. Dhabhar, N.E. Adler, J.D. Morrow, and
R.M. Cawthon. 2004. Accelerated telomere shortening in response to life stress.
\emph{Proceedings of the National Academy of Sciences (USA)} 101: 17312-17315.
}

\references{
\url{http://www.pnas.org/content/101/49/17312}
}

\examples{
data(Telomeres)
xyplot(years ~ telomere.length, Telomeres,
  xlab = "Time since diagnosis (years)",
  ylab = "Telomere length (ratio)"
)

plot(telomere.length ~ years, data = Telomeres,
  col = "red",
  pch = 16,
  xlab = "Chronicity (years)",
  ylab = "Telomere length (ratio)"
)
}

\keyword{datasets}
