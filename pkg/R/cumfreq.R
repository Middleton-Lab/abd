
panel.cumfreq <- function(x, type=c('smooth','step'), groups=NULL, ...) {
    if (!is.null(groups)) {
        panel.superpose(x, 
            ref = FALSE, groups = groups, 
            panel.groups = panel.cumfreq, 
            type = type, ...)
    }
    else {
		type=match.arg(type)
		if (type == 'step') {
			n <- length(x)
			xs <- rep(sort(x),each=2)
			p <- rep(ppoints(n-1),each=2)
			xs <- c(-Inf, xs, Inf)
			p <- c(0,0,p,1,1)
		} else {
			n <- length(x)
			xs <- sort(x)
			p <- ppoints(n)
			p <- aggregate(p, by=list(xs), mean)$x
			xs <- unique(xs)
		}
		panel.lines(x=xs, y=p, ...)
    }
}

cumfreq <- function(x, data, ...) { UseMethod('cumfreq') }

cumfreq.formula <- function(x, data=NULL, subscripts,
	...) {
	require(lattice)
	densityplot( x, data=data,
		ylab='cumulative frequency',
		panel=panel.cumfreq,
		prepanel=prepanel.cumfreq,
		...)
}

cumfreq.default <- function(x, ...) {
	cumfreq.formula( ~ x, ...)
}

prepanel.cumfreq <- function(x, ...) {
	list( xlim=range(x), ylim=c(0,1), dx=1, dy=1 )
}

