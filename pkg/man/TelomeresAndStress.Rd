\name{TelomeresAndStress}
\alias{TelomeresAndStress}
\docType{data}

\title{
Telomere Shortening
}

\description{
Telomere length (ratio) and years since their child's diagnosis with chronic illness.
}

\usage{data(TelomeresAndStress)}
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
Epel, E.S., E.H. Blackburn, J. Lin, F.S. Dhabhar, N.E. Adler, J.D. Morrow, and R.M. Cawthon. 2004. Accelerated telomere shortening in response to life stress. \emph{Proceedings of the National Academy of Sciences (USA)} 101: 17312-17315.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(TelomeresAndStress)

plot(telomere.length ~ years, data = TelomeresAndStress,
  col = "red",
  pch = 16,
  xlab = "Chronicity (years)",
  ylab = "Telomere length (ratio)")
}

\keyword{datasets}
