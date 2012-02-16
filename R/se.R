#' Standard Error of the Mean
#' 
#' \code{se} calculates the standard error of the mean of a vector.
#' 
#' 
#' @param x a vector of numeric values
#' @return A single numeric value of the standard error of the mean
#' @author Kevin M. Middleton (\email{kmm@@csusb.edu})
#' @seealso \code{\link{mean}} for the mean of a vector \code{\link{sd}} for
#' the standard deviation of a vector
#' @keywords univar
#' @export
#' @examples
#' 
#' set.seed(2)
#' n <- 10
#' y <- rnorm(n)
#' 
#' sd(y)
#' sd(y)/sqrt(n)
#' 
#' se(y)
#' 
se <- function(x){
  n <- length(x)
  sd(x)/sqrt(n)
  }
