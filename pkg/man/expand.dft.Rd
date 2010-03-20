\name{expand.dft}
\alias{expand.dft}
\title{
Expand a data.frame
}

\description{
This function expands a \code{data.frame} based on it's \code{Frequency} column. Useful for generating contingency tables from "summary" tables.
}

\usage{
expand.dft(x, na.strings = "NA", as.is = FALSE, dec = ".")
}

\arguments{
  \item{x}{a \code{data.frame}}
  \item{na.strings}{passed to \code{type.convert}}
  \item{as.is}{passed to \code{type.convert}}
  \item{dec}{passed to \code{type.convert}}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
  \item{DF}{a \code{data.frame} of \code{x} in "long" format.}
}

\references{
See https://stat.ethz.ch/pipermail/r-help/2009-January/185561.html for discussion of expand.dft().}

\author{
Modified for the \code{abd} package by Kevin Middleton from a function posted by Marc Schwartz to the r-help mailing list.
}

%\note{
%%  ~~further notes~~
%}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
data(AspirinCancer)
AspirinCancer

AspirinCancer.expanded <- expand.dft(AspirinCancer)
str(AspirinCancer.expanded)

# Plot 2X2 Contingency tables
plot( ~ Aspirin.treatment + Cancer, data = AspirinCancer.expanded)
plot(table(AspirinCancer.expanded), main = "")
}

\keyword{ manip }
