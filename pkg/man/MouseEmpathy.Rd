\name{MouseEmpathy}
\alias{MouseEmpathy}
\docType{data}

\title{
Mouse Empathy
}

\description{
Percentage of time spent stretching in three treatments of mice. Both \code{condition} and \code{treatment} code for the same variable.
}

\usage{data(MouseEmpathy)}
\format{
  A data frame with 42 observations on the following 3 variables.
  \describe{
    \item{\code{condition}}{a factor with levels \code{BothWrithing}, \code{Isolated}, and \code{OneWrithing}}
    \item{\code{percent.stretching}}{a numeric vector}
    \item{\code{treatment}}{a factor with levels \code{bw}, \code{isolated}, and \code{ow}}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Langford, D.J., S.E. Crager, Z. Shehzah, S.B. Smith, S.G. Sotocinal, J.S. Levenstadt, M.L. Chande, D.J. Levitin, J.S. Mogill. 2006. Social modulation of pain as evidence for empathy in mice. \emph{Science} 312: 1967-1970.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(MouseEmpathy)
str(MouseEmpathy)

aov.fit <- aov(percent.stretching ~ treatment, data = MouseEmpathy)
summary(aov.fit)
}

\keyword{datasets}
