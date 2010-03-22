\name{IndianRopeTrick}
\alias{IndianRopeTrick}
\docType{data}

\title{
Indian Rope Trick
}

\description{
Perceived impressiveness (\code{impressiveness.score}) of a written account of the Indian Rope Trick and the corresponding number of \code{year}s since it was witnessed.
}

\usage{data(IndianRopeTrick)}
\format{
  A data frame with 21 observations on the following 2 variables.
  \describe{
    \item{\code{years}}{a numeric vector}
    \item{\code{impressiveness.score}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Wiseman, R. and P. Lamont. 1996. Unravelling the Indian rope-trick. \emph{Nature} 383: 212-213.
}

\references{
\url{http://www.richardwiseman.com/resources/ropeJSPR.pdf}
}

\examples{
data(IndianRopeTrick)
IndianRopeTrick

rank.years <- c(1, 3.5, 3.5, 2, 5.5, 5.5, 13, 7, 8, 9, 10.5, 12, 14.5,
  17, 18, 19, 14.5, 10.5, 16, 20.5, 20.5)
rank.imp <- c(2, 2, 2, 5, 5, 5, 7, rep(12.5, times = 10),
  rep(19.5, times = 4))

sum.prods <- sum_of_products(rank.years, rank.imp)
SS.years <- sum_of_squares(rank.years)
SS.imp <- sum_of_squares(rank.imp)
sum.prods / (sqrt(SS.years) * sqrt(SS.imp))

# With cor.test(); Note warning about ties. See discussion on
# p. 446.
cor.test(IndianRopeTrick$years,
  IndianRopeTrick$impressiveness.score, method = "spearman")
}

\keyword{datasets}
