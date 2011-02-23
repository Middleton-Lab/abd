\name{Blackbirds}
\alias{Blackbirds}
\docType{data}

\title{
Testosterone Levels in Blackbirds
}

\description{
Experimental manipulation of testosterone levels in male Red-winged Blackbirds (\emph{Agelaius phoeniceus}) and resulting changes in antibody levels
}

\usage{data(Blackbirds)}

\format{
  A data frame with 13 observations on the following 6 variables.
  \describe{
    \item{\code{before}}{a numeric vector}
    \item{\code{after}}{a numeric vector}
    \item{\code{log.before}}{a numeric vector}
    \item{\code{log.after}}{a numeric vector}
    \item{\code{diff.in.logs}}{a numeric vector}
    \item{\code{diff}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Hasselquist, D., J.A. Marsh, P.W. Sherman, and J.C. Wingfield. 1999. 
Is avian immunocompetence suppressed by testosterone? \emph{Behavioral 
Ecology and Sociobiology} 45: 167-175.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(Blackbirds)
Blackbirds
xyplot(log.after ~ log.before, data = Blackbirds,
  ylab = "log Antibody production after implant",
  xlab = "log Antibody production before implant"
)
}

\keyword{datasets}
