repeatability <- function(vc){
  var.among <- vc[["var.among"]]
  var.within <- vc[["var.within"]] 
  R <- var.among / (var.among + var.within)
  class(R) <- "repeatability"
  R
}