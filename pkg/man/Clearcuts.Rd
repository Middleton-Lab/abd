\name{Clearcuts}
\alias{Clearcuts}
\docType{data}

\title{
Biomass Change in Rainforests near Clearcuts
}

\description{
Biomass change in 36 Amazonian rainforests following clearcuts ranging from 50m to several kilometers.
}

\usage{data(Clearcuts)}
\format{
  The format is:
 num [1:36] -10.8 -4.9 -2.6 -1.6 -3 -6.2 -6.5 -9.2 -3.6 -1.8 ...
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Laurance, W.F., S.G. Laurance, L.V. Ferreira, J.M. Rankin-de Merona,   C. Gascon, T.E. Lovejoy. 1997. Biomass collapse in Amazonian forest fragments. \emph{Science} 278: 1117-1118.
}

\references{
\url{http://www.sciencemag.org/cgi/content/abstract/278/5340/1117}
}

\examples{
data(Clearcuts)
str(Clearcuts)
hist(Clearcuts)
}

\keyword{datasets}
