# Data sets from *The Analysis of Biological Data* #

## Description ##

The `abd` package contains data sets and sample code for the book,
*The Analysis of Biological Data* by Michael C. Whitlock and Dolph
Schluter (2009; Roberts and Company Publishers).
   
## Authors ##

Kevin M. Middleton (CSU, San Bernardino); Randall Pruim (Calvin College)

## Reference ##

Whitlock, M.C. and D. Schluter. 2009. [*The Analysis of Biological 
Data*][abdurl]. Roberts and Company Publishers. ISBN: 0981519407.

[abdurl]:http://www.roberts-publishers.com/whitlock/

## Examples ##

     trellis.par.set(theme=col.abd())  # set color theme
     show.settings()
     findData(3)                       # look for data sets in chapter 3
     findData('Finch')                 # look for data sets with 'finch' in name
