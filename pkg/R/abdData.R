
abdData <- function(..., 
	chapters=1:21, 
	types= c('Example','Problem'),
	numbers=1:100, 
	pattern = '*',
	ignore.case=TRUE) {

	dots <- list(...)

	for (x in dots) {
		if ( all( x %in% c('Example','Problem') ) ) { 
			types <- x 
		} else { 
			if (is.character(x)) { 
				pattern <- x 
			} else { 
				if (is.numeric(x)) { 
					chapters <- x 
				} 
			}
		}
	}


	results <- subset(dataInfo,
				chapter %in% chapters & 
				type %in% types &
				number %in% numbers &
				grepl(pattern,name,ignore.case=ignore.case) 
	)

	if (prod(dim(results)) == 0) { return (NULL) }
	return(results)
}

