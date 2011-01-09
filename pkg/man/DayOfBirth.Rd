\name{DayOfBirth}
\alias{DayOfBirth}
\docType{data}

\title{
Day of Birth
}

\description{
Day of the week for 350 U.S. births in 1999.
}

\usage{data(DayOfBirth)}
\format{
  A data frame with 7 observations on the following 2 variables.
  \describe{
    \item{\code{day}}{a character vector}
    \item{\code{births}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Ventura, S.J., J.A. Martin, S.C. Curtin, F. Menacker, and B.E. Hamilton. 2001. Births: final data for 1999. \emph{National Vital Statistics Reports} Vol. 49, No. 1.
}

\references{
\url{http://cdc.gov/NCHS/products/nvsr.htm}
}

\examples{
data(DayOfBirth)
DayOfBirth
barchart( day ~ births, DayOfBirth)

# fix bad ordering of days
DayOfBirth$oday <- with(DayOfBirth, ordered(day, levels = day))
barchart( oday ~ births, DayOfBirth)
barchart( births ~ oday, DayOfBirth, horizontal = FALSE)
barchart( births ~ oday, DayOfBirth, horizontal = FALSE, 
 scales = list(x=list(rot=45)))

barplot(DayOfBirth$births,
  ylim = c(0, 70),
  names.arg = DayOfBirth$day,
  las = 2,
  mgp = c(3, 0.75, 0))

}

\keyword{datasets}
