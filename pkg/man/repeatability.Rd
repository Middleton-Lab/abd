\name{repeatability}
\alias{repeatability}

\title{
Repeatability
}

\description{
Calculate the repeatability of observations in a random effects ANOVA.
}

\usage{
repeatability(vc)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vc}{an object of class "varcomp"}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Returns a numeric value of the repeatability with a class "vc."
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
\code{\link{varcomps}}
}

\examples{
data(WalkingStickFemurs)
WalkingStickFemurs$specimen <- factor(WalkingStickFemurs$specimen)
aovfit <- aov(femurlength ~ 1 + Error(specimen), data = WalkingStickFemurs)
vc <- varcomps(aovfit, n = 2)
vc
repeatability(vc)
}

\keyword{ univar }

