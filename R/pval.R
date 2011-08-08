# Contributed by Randall Pruim (rpruim@calvin.edu)

pval <- function(x, ...){UseMethod("pval", x)}

pval.htest <- function (x, digits=4,...){
  pval <- x$p.value
  stat <- x$statistic
  param <- x$parameter
  alt <- x$alternative
  method <- x$method
  null <- x$null.value
  estimate <- x$estimate
  direction <- switch(alt, 
  	'less' = ' < ',
  	'greater' = ' > ',
  	'two.sided' = ' <> '
	)
  cat('\n')
  cat(paste('Method: ', method,  sep=""))
  cat('\n\n')
  cat(paste(
  	'Null Hypothesis: ', 
  	names(null), 
	" = ", 
	null,
  	sep="") 
  )  
  cat('\n')
  cat(paste(
  	'Alt. Hypothesis: ', 
  	names(null), 
	direction, 
	null,
  	sep="") 
  )  
  cat('\n\n')
  cat(paste(names(stat), " = ", 
  	signif(stat,digits=digits),
  	sep="") )  
  cat('  (')
  cat( paste( 
  	names(param), " = ", 
  	signif(param,digits=digits), 
  	sep="",
	collapse=', ') )  
  cat(')\n\n')
  cat( paste("p-value = ", signif(pval,digits), sep="") ) 
  cat('\n\n')
  invisible(pval)
}
