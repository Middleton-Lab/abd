# see https://stat.ethz.ch/pipermail/r-help/2009-January/185561.html
# for discussion of expand.dft(). Modified for abd package.

#V3
expand.dft <- function(x, col.exp, na.strings = "NA", as.is = FALSE, dec = ".") {
  # Error checking
  if (missing(x)) {
    stop("x is missing. You must supply a data.frame.")
  }
  if (missing(col.exp)) {
    stop("col.exp is missing You must supply a column name to expand.")
  }
  if (class(x) != "data.frame") {
    stop("x must be a data.frame.")
  }
  if (class(col.exp) != "character") {
  	  if (class(col.exp) != "numeric") {
        stop("col.exp must be character or numeric.")
      }
   }
  
  # Get the column name if col.exp is numeric
  if (class(col.exp) == "numeric") col.exp <- names(x)[col.exp]
  
  # Frequencies to expand via rep()
  to.expand <- x[, col.exp]
  
  DF <- sapply(1:nrow(x), function(i) x[rep(i, each = to.expand[i]), ], 
    simplify = FALSE)
  DF <- subset(do.call("rbind", DF), select = -(get(col.exp)))
  for (i in 1:ncol(DF)) {
    DF[[i]] <- type.convert(as.character(DF[[i]]),
      na.strings = na.strings, as.is = as.is, dec = dec)
  }
DF
} 

#V2
#expand.dft <- function(x, col.exp, na.strings = "NA", as.is = FALSE, dec = "."){
#  to.expand <- eval(substitute(x[a], list(a = col.exp)))[, 1]
#  DF <- sapply(1:nrow(x), function(i) x[rep(i, each = to.expand[i]), ], 
#    simplify = FALSE)
#  DF <- subset(do.call("rbind", DF), select = -(get(col.exp)))
#  for (i in 1:ncol(DF)){
#    DF[[i]] <- type.convert(as.character(DF[[i]]),
#      na.strings = na.strings, as.is = as.is, dec = dec)
#    }
#DF
#} 

# V1
#expand.dft <- function(x, col.exp, na.strings = "NA", as.is = FALSE, dec = "."){
#  
#  DF <- sapply(1:nrow(x), function(i) x[rep(i, each = x$Frequency[i]), ], 
#    simplify = FALSE)
#  DF <- subset(do.call("rbind", DF), select = -Frequency)
#  for (i in 1:ncol(DF)){
#    DF[[i]] <- type.convert(as.character(DF[[i]]),
#      na.strings = na.strings, as.is = as.is, dec = dec)
#    }
#DF
#} 
