varcomps <- function(fm, n){
  summ <- summary(fm)
  MS.groups <- summ[[1]][[1]]$"Mean Sq"
  MS.error <- summ[[2]][[1]]$"Mean Sq"
  F <- MS.groups / MS.error
  df1 <- summ[[1]][[1]]$"Df"
  df2 <- summ[[2]][[1]]$"Df"
  p <- pf(F, df1, df2, lower.tail = FALSE)
  
  var.within <- MS.error
  var.among <- (MS.groups - MS.error) / n
  
  varcomp.obj <- list("MS.groups" = MS.groups,
                      "MS.error" = MS.error,
                      "F" = F,
                      "df1" = df1,
                      "df2" = df2,
                      "p" = p,
                      "var.within" = var.within,
                      "var.among" = var.among)
  
  class(varcomp.obj) <- "varcomps"
  varcomp.obj
}