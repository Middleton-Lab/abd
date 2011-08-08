sum_of_products <- function(x, y){
  sum((x - mean(x)) * (y - mean(y)))
}

sum_of_squares <- function(x){
  sum((x - mean(x))^2)
}