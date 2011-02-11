varcomps <- function(fm, n){
  summ <- summary(fm)
  MS.groups <- summ[[1]][[1]]$"Mean Sq"
  MS.error <- summ[[2]][[1]]$"Mean Sq"
  Fstat <- MS.groups / MS.error
  df1 <- summ[[1]][[1]]$"Df"
  df2 <- summ[[2]][[1]]$"Df"
  p <- pf(Fstat, df1, df2, lower.tail = FALSE)
  
  var.within <- MS.error
  var.among <- (MS.groups - MS.error) / n
  
  varcomp.obj <- list("MS.groups" = MS.groups,
                      "MS.error" = MS.error,
                      "Fstat" = Fstat,
                      "df1" = df1,
                      "df2" = df2,
                      "p" = p,
                      "var.within" = var.within,
                      "var.among" = var.among)
  class(varcomp.obj) <- "varcomps"
  return(varcomp.obj)
}

print.varcomps <- function(x, digits = 3, ...){
  cat("Mean Squares\n")
  cat("\tGroups\t", x$MS.groups, "\n")
  cat("\tError\t", x$MS.error, "\n")
  cat("Variance Components\n")
  cat("\tWithin\t", x$var.within, "\n")
  cat("\tAmong\t", x$var.among, "\n")
  cat("F = ", signif(x$Fstat, digits = digits), "\n")
  cat("p = ", signif(x$p, digits = digits), "on", x$df1, "and", x$df2, "degrees of freedom.\n")
}
