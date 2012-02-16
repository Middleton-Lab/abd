#' Sum of Squares and Sum of Products
#' 
#' Calculate the sum of squares (\code{sum_of_squares}) or sum of products
#' (\code{sum_of_products}) of a numeric vector or pair of numeric vectors
#' 
#' 
#' @rdname sum-of-products
#' @param x a vector of numeric values
#' @param y a vector of numeric values
#' @return a numeric value
#' @author Kevin M. Middleton (\email{kmm@@csusb.edu})
#' @keywords univar
#' @export
#' @examples
#' 
#' set.seed(4)
#' x <- rnorm(10)
#' sum_of_squares(x)
#' 
#' y <- rnorm(10)
#' sum_of_products(x, y)
#' 
sum_of_squares <- function(x){
  sum((x - mean(x))^2)
}

#' @rdname sum-of-products
#' @export
sum_of_products <- function(x, y){
  sum((x - mean(x)) * (y - mean(y)))
}



