\name{HornedLizards}
\alias{HornedLizards}

\docType{data}

\title{
Horn Length and Predation Status of Horned Lizards
}
\description{
Squamosal horn length (mm; \code{horn.length}) and predation status
(\code{group}; living or killed) for 184 horned lizards
(\emph{Phrynosoma mcalli}).
}

\usage{data(HornedLizards)}

\format{
  A data frame with 184 observations on the following 2 variables.
  \describe{
    \item{\code{horn.length}}{a numeric vector}
    \item{\code{group}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Young, K.V., E.D. Brodie, Jr., and E.D. Brodie, III. 2004. How the horned lizard got its horns. \emph{Science} 304: 65.
}

\references{
\url{http://www.sciencemag.org/cgi/pdf_extract/304/5667/65}
}

\examples{
data(HornedLizards)
str(HornedLizards)

histogram(~horn.length | group, HornedLizards, 
	layout=c(1,2),
	xlab="Horn Length (mm)")
}

\keyword{datasets}
