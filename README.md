<!-- badges: start -->
[![RStudio_CRAN_mirror_downloads_badge](http://cranlogs.r-pkg.org/badges/grand-total/abd?color=blue)](https://cran.r-project.org/web/packages/abd)
<!-- badges: end -->
  
# Data sets from The Analysis of Biological Data

## Description

The `abd` package contains data sets and sample code for the book,
*The Analysis of Biological Data* by Michael C. Whitlock and Dolph
Schluter (2009; Roberts and Company Publishers).
   
## Authors

Kevin M. Middleton (University of Missouri); Randall Pruim (Calvin College)

## Reference

Whitlock, M.C. and D. Schluter. 2009. [*The Analysis of Biological 
Data*][abdurl]. Roberts and Company Publishers. ISBN: 0981519407.

[abdurl]:https://whitlockschluter.zoology.ubc.ca/

## Installation

You can install the most recent version directly from github using
`install_github()` from the [remotes
package](https://cran.r-project.org/package=remotes).

```R
remotes::install_github("Middleton-Lab/abd")
```

Note that if you are using Windows, you may first need to install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/.

Stable release versions of `abd` are available via CRAN:

```R
install.packages("abd")
```

## Examples

```R
trellis.par.set(theme=col.abd())  # set color theme
show.settings()
abdData(3)                        # look for data sets in chapter 3
abdData('Finch')                  # look for data sets with 'finch' in name
```
