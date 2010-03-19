\name{varcomps}
\alias{varcomps}

\title{
Variance Components
}

\description{
Extract variance components from a mixed-model ANOVA.
}

\usage{
varcomps(fm, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fm}{an object of class "aov"}
  \item{n}{an integer giving the within-group sample size}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
List of class "varcomps" with the following components:
\item{MS.groups}{xxx}
\item{MS.error}{xxx}
\item{F}{xxx}
\item{df1}{xxx}
\item{df2}{xxx}
\item{p}{xxx}
\item{var.within}{xxx}
\item{var.among}{xxx}
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
\code{\link{repeatability}}
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
