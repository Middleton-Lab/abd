\name{expand.dft}
\alias{expand.dft}
\title{
Expand a data.frame
}

\description{
This function expands a \code{data.frame} based on a column (\code{col.exp}). Useful for generating long tables from contingency tables.
}

\usage{
  expand.dft(x, col.exp, na.strings = "NA", as.is = FALSE, dec = ".")
}

\arguments{
  \item{x}{a \code{data.frame}}
  \item{col.exp}{a character string of the column name or the column number to use for "expansion"}
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
See \url{https://stat.ethz.ch/pipermail/r-help/2009-January/185561.html} for discussion of \code{expand.dft()}.}

\author{
Modified for the \code{abd} package by Kevin Middleton from a function posted by Marc Schwartz to the r-help mailing list.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
  \code{\link{table}}
}

\examples{
data(AspirinCancer)
AspirinCancer

# Specifying col.exp as character
AspirinCancer.expanded <- expand.dft(AspirinCancer, "Frequency")
str(AspirinCancer.expanded)

# Specifying col.exp as numeric
AspirinCancer.expanded <- expand.dft(AspirinCancer, 3)
str(AspirinCancer.expanded)

# Plot 2X2 Contingency tables
plot( ~ Aspirin.treatment + Cancer, data = AspirinCancer.expanded)
plot(table(AspirinCancer.expanded), main = "")
}

\keyword{ manip }
