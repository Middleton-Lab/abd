
sample <- function (x, size, replace=FALSE, prob=NULL) { UseMethod('sample') }

sample.default <- function(x, ...) { base::sample(x, ...) }

sample.data.frame <- function(x, size, replace = FALSE, prob = NULL) {
	n <- nrow(x)
	ids <- base::sample(n, size, replace=replace, prob=prob)
	data <- data.frame( x [ ids, ] )
	names(data) <- names(x)
	data$orig.row <- ids
	if (size < 50) { return(data) } else {invisible(data)}
}
