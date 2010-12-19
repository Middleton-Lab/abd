\name{interval}
\alias{interval}
\alias{interval.htest}

\title{Confidence Interval}
\description{
Extract a confidence interval from an htest object.
}

\usage{interval(x, ...)}

\arguments{
  \item{x}{An object of class \code{htest}.}
  \item{\dots}{Additional arguments.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Randall Pruim (rpruim@calvin.edu)
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
interval(t.test(rnorm(100)))
}

\keyword{ univar }
