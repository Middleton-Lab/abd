\name{LizardBite}
\alias{LizardBite}
\docType{data}

\title{
Bite Force in Collard Lizards
}

\description{
Bite force (N) and territory area in 11 male collared lizards (\emph{Crotaphytus collaris}).
}

\usage{data(LizardBite)}
\format{
  A data frame with 11 observations on the following 2 variables.
  \describe{
    \item{\code{bite}}{force of bite (N)}
    \item{\code{territory}}{area of territory}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Lappin, A. K., and J. F. Husak. 2005. Weapon performance, not size, determines mating success and potential reproductive output in the collared lizard (\emph{Crotaphytus collaris}). \emph{The American Naturalist} 166: 426-436.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(LizardBite)
str(LizardBite)
xyplot(territory ~ bite, LizardBite)
}

\keyword{datasets}
