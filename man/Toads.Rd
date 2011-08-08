\name{Toads}
\alias{Toads}
\docType{data}

\title{
Right-handed Toads
}

\description{
Hypothetical probability of a toad being right-handed
}

\usage{data(Toads)}
\format{
  A data frame with 19 observations on the following 2 variables.
  \describe{
    \item{\code{n.toads}}{a numeric vector}
    \item{\code{prob}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\examples{
data(Toads)
Toads
xyplot(prob~n.toads, Toads, type='h', lwd=4)
barchart(prob~n.toads, Toads, horizontal=FALSE)

}

\keyword{datasets}
