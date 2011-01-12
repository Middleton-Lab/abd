propCI <- function(..., 
	method=c('binom.test','prop.test','wilson','agresti-coull')) 
{
	if (is.character(method)) {
		method <- match.arg(method)
	}
	switch( method,
		prop.test = interval(prop.test(...)),
		binom.test = interval(binom.test(...)),
		wilson = wilsonCI(...),
		'agresti-coull' = wilsonCI(...)
	)
}
