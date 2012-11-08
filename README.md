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

[abdurl]:http://www.roberts-publishers.com/whitlock/

## Installation

You can install the most recent version directly from github using
`install_github()` from the [devtools
package](https://github.com/hadley/devtools).

```R
require(devtools)
install_github("abd", "kmiddleton")
```

Note that if you are using Windows, you will first need to install
Rtools. Start at <http://www.murdoch-sutherland.com/Rtools/> and
follow the links to CRAN.

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
