\name{repeatability}

\alias{repeatability}
\alias{print.repeatability}

\title{
Repeatability
}

\description{ 
Calculate the repeatability (intraclass correlation coefficient) of observations in a random effects ANOVA fit by \code{\link{aov}} with \code{Error()} via the \code{\link{varcomps}} functon, a one-way ANOVA fit with \code{aov()}, or a or a linear mixed effects model fit by \code{\link{lme}}. The first two methods follow Lessels and Boag (1987).
}

\usage{
repeatability(x)
}

\arguments{
  \item{x}{an object of class \code{varcomps} (from \code{\link{varcomps}}),
  class \code{aov},  or \code{lme} (from \code{\link{lme}})}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
	\code{repeatability} returns a numeric value of the repeatability (intraclass correlations coefficient) with a \link{class} \code{repeatability} and \code{rep.aov} \emph{or} \code{rep.lme}, depending on the 	class of \code{x}.
}

\references{
Lessels C.M, and P.T. Boag. 1987. Unrepeatable repeatabilities: a common mistake. 
\emph{Auk} 104: 116-121.
}

\author{
Kevin M. Middleton (\email{kmm@csusb.edu})
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
# With aov() and Error()
Error.fit <- aov(femur.length ~ 1 + Error(specimen), data = WalkingStickFemurs)
vc <- varcomps(Error.fit, n = 2)
vc
repeatability(vc)

# With aov()
aov.fit <- aov(femur.length ~ specimen, data = WalkingStickFemurs)
repeatability(aov.fit)

# With lme()
lme.fit <- lme(femur.length ~ 1, random = ~ 1 | specimen, 
               data = WalkingStickFemurs)
repeatability(lme.fit)
}

\keyword{ univar }

