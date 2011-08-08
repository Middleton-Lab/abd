propPval <- function(..., method=c('binom.test', 'prop.test')) {
	if (is.character(method)) {
		method <- match.arg(method)
	}
	switch( method,
		prop.test = pval(prop.test(...)),
		binom.test = pval(binom.test(...)),
	)
}
