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
data(Aspirin)
Aspirin

# Specifying col.exp as character
Aspirin.expanded <- expand.dft(Aspirin, "count")
str(Aspirin.expanded)
xtabs(~treatment + cancer, Aspirin.expanded)

# Specifying col.exp as numeric
Aspirin.expanded <- expand.dft(Aspirin, 3)
str(Aspirin.expanded)
xtabs(~treatment + cancer, Aspirin.expanded)

# Plot 2X2 Contingency tables
plot( ~ treatment + cancer, data = Aspirin.expanded)
plot(table(Aspirin.expanded), main = "")
mosaicplot(~treatment + cancer, Aspirin.expanded)

# much nicer looking plots using vcd
if(require(vcd)) {
	mosaic(~treatment + cancer, Aspirin.expanded)
}
}

\keyword{ manip }
