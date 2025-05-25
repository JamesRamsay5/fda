
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fda

<!-- badges: start -->

[![R-CMD-check](https://github.com/JamesRamsay5/fda/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JamesRamsay5/fda/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/fda)](https://CRAN.R-project.org/package=fda)
<!-- badges: end -->

The goal of fda is to …

## Installation

You can install the development version of fda from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JamesRamsay5/fda")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fda)
#> Loading required package: splines
#> Loading required package: fds
#> Loading required package: rainbow
#> Loading required package: MASS
#> Loading required package: pcaPP
#> Loading required package: RCurl
#> Loading required package: deSolve
#> 
#> Attaching package: 'fda'
#> The following object is masked from 'package:graphics':
#> 
#>     matplot
#> The following object is masked from 'package:datasets':
#> 
#>     gait
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
