\name{Powerball}
\alias{Powerball}
\docType{data}

\title{
Powerball Tickets Sold
}

\description{
The number of Powerball tickets sold per day of the week for three years.
}

\usage{data(Powerball)}
\format{
  A data frame with 7 observations on the following 2 variables.
  \describe{
    \item{\code{day}}{a character vector}
    \item{\code{millions.of.tickets.sold}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Oster, E. 2004. Dreaming big: Why do people play Powerball? \emph{Chance News} 13.02.
}

\references{
\url{http://www.dartmouth.edu/~chance/chance_news/recent_news/chance_news_13.02.html}
}

\examples{
data(Powerball)
Powerball
xyplot(millions.of.tickets.sold ~ day, Powerball)
}

\keyword{datasets}
