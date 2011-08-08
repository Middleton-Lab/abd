favstats <-
function (x, na.rm = TRUE) 
{
    val <- c(
		median(x, na.rm=na.rm),
		IQR(x, na.rm=na.rm),
		mean(x, na.rm = na.rm), 
		sd(x, na.rm = na.rm), 
        var(x, na.rm = na.rm)
		)
    names(val) <- c("median", "IQR", "mean", "sd", "var")
    return(val)
}
