\docType{data}
\name{Powerball}
\alias{Powerball}
\title{Powerball Tickets Sold}
\format{A data frame with 7 observations on the following 2 variables.
\describe{ \item{day}{a character vector}
\item{millions.of.tickets.sold}{a numeric vector} }}
\source{
  Oster, E. 2004. Dreaming big: Why do people play
  Powerball? \emph{Chance News} 13.02.
}
\description{
  The number of Powerball tickets sold per day of the week
  for three years.
}
\examples{
data(Powerball)
Powerball
xyplot(millions.of.tickets.sold ~ day, Powerball)
}
\references{
  \url{http://www.dartmouth.edu/~chance/chance_news/recent_news/chance_news_13.02.html}
}
\keyword{datasets}

