sum.of.products <- function(x, y){
  sum((x - mean(x)) * (y - mean(y)))
}

sum.of.squares <- function(x){
  sum((x - mean(x))^2)
}