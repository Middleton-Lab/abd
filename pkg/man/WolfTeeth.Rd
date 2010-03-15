\name{WolfTeeth}
\alias{WolfTeeth}
\docType{data}
\title{
Wolf Tooth Measurements
}

\description{
Measurement (cm) of the distance between the canine and last molar teeth in 35 wolves.
}

\usage{data(WolfTeeth)}

\format{
  The format is:
 num [1:35] 10.2 10.4 9.9 10.7 10.3 9.7 10.3 10.7 10.1 10.6 ...
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Whitlock, M. 1996. The heritability of fluctuating asymmetry and the genetic control of developmental stability. \emph{Proceedings of the Royal Society, Series B} 263: 849-853.
}

\references{
\url{http://rspb.royalsocietypublishing.org/content/263/1372/849.abstract}
}

\examples{
data(WolfTeeth)
WolfTeeth
hist(WolfTeeth)
}

\keyword{datasets}
