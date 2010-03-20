\name{repeatability}

\alias{repeatability}
\alias{print.repeatability}

\title{
Repeatability
}

\description{
Calculate the repeatability of observations in a random effects ANOVA fit by \code{\link{aov}} or a or a linear mixed effects model fit by \code{\link{lme}}.
}

\usage{
repeatability(x)
}

\arguments{
  \item{x}{an object of class \code{varcomps} (from \code{\link{varcomps}}) or \code{lme} (from \code{\link{lme}})}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\code{repeatability} returns a numeric value of the repeatability with a \link{class} \code{repeatability} and \code{rep.aov} \emph{or} \code{rep.lme}, depending on the class of \code{vc}.
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Kevin M. Middleton (kmm@csusb.edu)
}

%\note{
%%  ~~further notes~~
%}

\seealso{
  \code{\link{varcomps}},
  \code{\link{lme}}
}

\examples{
data(WalkingStickFemurs)
aovfit <- aov(femurlength ~ 1 + Error(specimen), data = WalkingStickFemurs)
vc <- varcomps(aovfit, n = 2)
vc
R.varcomps <- repeatability(vc)
R.varcomps
}

\keyword{ univar }

