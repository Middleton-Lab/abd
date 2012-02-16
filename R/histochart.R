#' Histogram from tabulated data
#' 
#' Uses \code{\link{barchart}} to build a histogram from tabulated data.
#' 
#' This is just a convenience wrapper around \code{\link{barchart}}.
#' 
#' @param x formula of form \code{frequency ~ value}
#' @param data data frame in which the formula \code{x} is interpreted
#' @param box.ratio ratio of bar widths to gaps between bars
#' @param origin where do bars begin?
#' @param horizontal Should bars go horizontal?
#' @param \dots other arguments passed to \code{\link{barchart}}
#' @author Randall Pruim (\email{rpruim@@calvin.edu})
#' @seealso \code{\link{barchart}}
#' @keywords graphics
#' @export
#' @examples
#' 
#' histochart( dbinom(0:30, 30, 0.35) ~ 0:30 )
#' 
histochart <- function(x, data=NULL, box.ratio=100, origin=0,
	horizontal=FALSE, ...) 
{
	barchart(x, data, box.ratio=box.ratio, origin=origin,
		horizontal=horizontal, ...)
}
