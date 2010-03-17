\name{StalkieEyespan}
\alias{StalkieEyespan}
\docType{data}

\title{
Stalk-eyed Fly Eyespan
}

\description{
Eyespan width (mm; \code{Eye.span}) in 45 stalk-eyed flies (\emph{Cyrtodiopsis dalmanni}) fed a corn or cotton diet \code{Food}.
}

\usage{data(StalkieEyespan)}
\format{
  A data frame with 45 observations on the following 2 variables.
  \describe{
    \item{\code{Food}}{a factor with levels \code{Corn} \code{Cotton}}
    \item{\code{Eye.span}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
David, P., T. Bjorksten, K. Fowler, and A. Pomiankowski. 2000. Condition-dependent signalling of genetic variation in stalk-eyed flies. \emph{Nature} 406: 186-188.
}

\references{
\url{http://www.nature.com/nature/journal/v406/n6792/full/406186a0.html}
}

\examples{
data(StalkieEyespan)
StalkieEyespan
str(StalkieEyespan)
}

\keyword{datasets}
