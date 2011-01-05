\name{Mosquitoes}
\alias{Mosquitoes}
\docType{data}

\title{
Body Size in Anopheles Mosquitoes
}

\description{
Weights of female and male mosquitos (\emph{Anopheles darlingi})
}

\usage{data(Mosquitoes)}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{weight}}{a numeric vector}
    \item{\code{sex}}{a factor with levels \code{female} and \code{male}}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Lounibos, L.P., N. Nishimura, J. Conn, and R. Lourenco-de-Oliveira. 1995. Life history correlates of adult size in the malaria vector \emph{Anopheles darlingi}. \emph{Memórias do Instituto Oswaldo Cruz} 90: 769-774.
}

\references{
\url{http://www.bioline.org.br/request?oc95154}
}

\examples{
data(Mosquitoes)
xyplot(weight ~ sex, Mosquitoes)
}

\keyword{datasets}
