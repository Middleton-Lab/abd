\name{Trematodes}
\alias{Trematodes}
\docType{data}

\title{
Frequencies of Fish Eaten by Trematode Infection Level
}

\description{
Frequencies of killifish (\emph{Fundulus parvipinnis}) eaten by birds depending
on level of infection by the trematode \emph{Euhaplorchis californiensis}.
}

\usage{data(Trematodes)}
\format{
  A data frame with 141 observations on the following 2 variables.
  \describe{
    \item{\code{infection.status}}{a factor with levels: \code{high}, \code{light}, and \code{uninfected}}
    \item{\code{eaten}}{a factor with levels: \code{no} and \code{yes}}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Lafferty, K.D. and A.K. Morris. 1996. Altered behavior of parasitized killifish increases susceptibility to predation by bird final hosts. \emph{Ecology} 77: 1390-1397.
}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(Trematodes)
xtabs(~ infection.status + eaten, Trematodes)
chisq.test( xtabs(~ infection.status + eaten, Trematodes) )
summary(chisq.test( xtabs(~ infection.status + eaten, Trematodes) ) )
}


\keyword{datasets}
