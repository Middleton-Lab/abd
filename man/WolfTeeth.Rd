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
A data frame with 35 observations of one variable.  
\describe{
\item{\code{length}}{distance from canine to last molar teach (cm)}
}
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
histogram(~length, WolfTeeth)
}

\keyword{datasets}
