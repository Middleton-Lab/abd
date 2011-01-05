\name{varcomps}

\alias{varcomps}
\alias{print.varcomps}

\title{
Variance Components
}

\description{
Extract variance components from a mixed-model ANOVA (i.e., fit with \code{+ Error(...)}).
}

\usage{
varcomps(fm, n)
}

\arguments{
  \item{fm}{an fitted model of class \code{aov}}
  \item{n}{an integer giving the within-group sample size}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\code{varcomps} returns a list of \link{class} \code{varcomps} with the following components:
\item{MS.groups}{Group mean square}
\item{MS.error}{Error mean square}
\item{F}{\emph{F}-statistic}
\item{df1}{Numerator degrees of freedom}
\item{df2}{Denominator degrees of freedom}
\item{p}{\emph{p}-value}
\item{var.within}{Within-group variance; equals \code{MS.error}}
\item{var.among}{Among group variance}
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Kevin M. Middleton (\email{kmm@csusb.edu})
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{repeatability}}
}

\examples{
data(WalkingStickFemurs)
aovfit <- aov(femur.length ~ 1 + Error(specimen), data = WalkingStickFemurs)
vc <- varcomps(aovfit, n = 2)
vc
R.varcomps <- repeatability(vc)
R.varcomps
}

\keyword{ univar }
