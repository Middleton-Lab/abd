as.xtabs <- function(x, ...) { UseMethod('as.xtabs') }

as.xtabs.data.frame <- function(x, rowvar=NULL, colvar=NULL, labels=1, ...) {

  if (labels >= 1) {
  	cnames <- names(x)[-1]
  	m <- as.matrix(x[,-c(labels)])
  } else {
  	cnames <- names(x)
	m <- as.matrix(x)
  }

  rnames <- x[,labels]
  rownames(m) <- rnames
  if (! is.character(rowvar) ) { rowvar <- "variable.1" }
  if (! is.character(colvar) ) { colvar <- "variable.2" }
  
  dn <- list( rnames, cnames)
  names(dn) <- c(rowvar, colvar)
  attr(m,'dimnames') <- dn
  class(m) <- c('xtabs', 'table')
  return(m)
}

as.xtabs.matrix <- function(x, rowvar=NULL, colvar=NULL, ...)  {
  rnames <- rownames(x)
  cnames <- colnames(x)
  rownames(m) <- rnames
  if (! is.character(rowvar) ) { rowvar <- "variable.1" }
  if (! is.character(colvar) ) { colvar <- "variable.2" }
  
  dn <- list( rnames, cnames)
  names(dn) <- c(rowvar, colvar)
  attr(m,'dimnames') <- dn
  class(m) <- c('xtabs', 'table')
  return(m)
}
