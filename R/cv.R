#' Coefficient of Variation
#' 
#' Compute coefficient of variation (ratio of standard deviation to mean).
#' 
#' 
#' @param x A numeric vector
#' @author Randall Pruim (\email{rpruim@@calvin.edu})
#' @keywords univar stats
#' @export
#' @examples
#' 
#' cv(GlidingSnakes$undulation.rate)
#' 
cv <- function(x) {
	sd(x) / mean(x)
}
