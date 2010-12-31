
Sample <- function (x, size, replace=FALSE, prob=NULL) { UseMethod('Sample') }

Sample.default <- function(x, ...) { sample(x, ...) }

Sample.data.frame <- function(x, size, replace = FALSE, prob = NULL) {
	n <- nrow(x)
	ids <- sample(n, size, replace=replace, prob=prob)
	data <- data.frame( x [ ids, ] )
	names(data) <- names(x)
	data$orig.row <- ids
	if (size < 50) { return(data) } else {invisible(data)}
}
