[![CRAM_Status_Badge](http://www.r-pkg.org/badges/version/vcpen)](https://CRAN.R-project.org/package=smartDesign)
[![Downloads](http://cranlogs.r-pkg.org/badges/vcpen)](https://CRAN.R-project.org/package=smartDesign)
[![Total-Downloads](https://cranlogs.r-pkg.org/badges/grand-total/smartdes)](https://CRAN.R-project.org/package=smartDesign)

# The `smartDesign` Package
Sequential Multiple Assignment Randomized Trial Design, with two main design options:
Single Sequential Treatment (SST) and Dynamic Treatment Regminen (DTR) designs, for which we split into separate calculations and plot methods.

# SST

## `smartSST()`

Performs sensitivity and specificity calculations for the SST design given means for G0 and G1, and sample size. User specifies the B level to be applied in the power and the plot functions.

## `powerSST()`
Power for the SST design.

# DTR

## `smartDTR()`

Sensitivity and Specificity calculations for given mean and variance vectors for the Dynamic Treatment Regimen (DTR).

## `powerDTR()`
Power to compere two groups of DTR settings. Plot and print methods are available.

