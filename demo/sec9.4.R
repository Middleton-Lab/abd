VampireBites

## Table 9.4-1
xtabs(count ~ estrous + bitten, data = VampireBites)

## Fisher's Exact Test
fisher.test(xtabs(count ~ estrous + bitten, data = VampireBites))

## Section 9.5
## Using G-test
## Source from http://www.psych.ualberta.ca/~phurd/cruft/
try({
  source("http://www.psych.ualberta.ca/~phurd/cruft/g.test.r")
  print( g.test(xtabs(count ~ estrous + bitten, data = VampireBites)) )
  cat("\nObserved Values:\n")
  print(xtabs(count ~ estrous + bitten, data = VampireBites))
  cat("\nExpected Values:\n")
  g.test(xtabs(count ~ estrous + bitten, data = VampireBites))$expected
  }
)
