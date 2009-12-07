\name{WillsPresidents}
\alias{WillsPresidents}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(WillsPresidents)}
\format{
  A data frame with 16 observations on the following 5 variables.
  \describe{
    \item{\code{Candidate}}{a factor with levels \code{Bush} \code{Carter} \code{Clinton} \code{Dole} \code{Dukakis} \code{Ford} \code{Gore} \code{Kennedy} \code{Kerry} \code{Nixon} \code{Reagan}}
    \item{\code{Winner.}}{a factor with levels \code{n} \code{y}}
    \item{\code{Rounded.wills}}{a numeric vector}
    \item{\code{Loser}}{a factor with levels \code{.} \code{122} \code{149} \code{18} \code{32} \code{33} \code{48} \code{75} \code{85}}
    \item{\code{Difference}}{a factor with levels \code{.} \code{1} \code{20} \code{23} \code{26} \code{27} \code{36} \code{4} \code{41}}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(WillsPresidents)
## maybe str(WillsPresidents) ; plot(WillsPresidents) ...
}
\keyword{datasets}
