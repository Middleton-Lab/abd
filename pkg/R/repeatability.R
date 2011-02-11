repeatability <- function(x){
  if(inherits(x, "varcomps")){
    var.among <- x[["var.among"]]
    var.within <- x[["var.within"]] 
    R <- var.among / (var.among + var.within)
    class(R) <- c("repeatability", "rep.varcomps")
  }
  if(inherits(x, "aov")){
    smry <- summary(x)
    n <- smry[[1]]$Df[1] + 1
    MSw <- smry[[1]]$"Mean Sq"[2]
    MSa <- smry[[1]]$"Mean Sq"[1]
    s2a <- (MSa - MSw)/n					
    R <- s2a/(MSw + s2a)
    class(R) <- c("repeatability", "rep.aov")
  }
  if(inherits(x, "lme")){
    varcomps <- VarCorr(x)
    var.among <- as.numeric(varcomps[1, 1])
    var.within <- as.numeric(varcomps[2, 1])
    R <- var.among / (var.among + var.within)
    class(R) <- c("repeatability", "rep.lme")
  }
  return(R)
}

print.repeatability <- function(x, ...){
  if (inherits(x, "rep.varcomps")){
    cat("Repeatability measured by random effects ANOVA.\n")
    cat("\tRepeatability is", round(x, digits = 3), "\n")
  } 
  if (inherits(x, "rep.aov")){
    cat("Repeatability measured by one-way ANOVA.\n")
    cat("\tRepeatability is", round(x, digits = 3), "\n")
  } 
  if (inherits(x, "rep.lme")){
    cat("Repeatability measured by linear mixed-effects model.\n")
    cat("\tRepeatability is", round(x, digits = 3), "\n")
  }
}
