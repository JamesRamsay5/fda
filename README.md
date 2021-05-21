
<!-- README.md is generated from README.Rmd. Please edit then 'knit' that file -->

# fda

<!-- badges: start -->

[![R-CMD-check](https://github.com/JamesRamsay5/fda/workflows/R-CMD-check/badge.svg)](https://github.com/JamesRamsay5/fda/actions)
<!-- badges: end -->

The ‘fda’ package supports functional data analysis in R, as described
in [Ramsay, Hooker, and Graves (2009) Functional Data Analysis with R
and MATLAB
(Springer)](https://www.amazon.com/Functional-Data-Analysis-MATLAB-Use/dp/0387981845)

## Installation

You can install the released version of fda from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fda")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JamesRamsay5/fda")
```

## Example

The `fda` package comes with script files that run nearly all the
examples in Ramsay, Hooker and Graves (2009). You can find the script
files as follows:

``` r
library(fda)
#> Loading required package: splines
#> Loading required package: Matrix
#> Loading required package: fds
#> Loading required package: rainbow
#> Loading required package: MASS
#> Loading required package: pcaPP
#> Loading required package: RCurl
#> 
#> Attaching package: 'fda'
#> The following object is masked from 'package:graphics':
#> 
#>     matplot
(scripts <- system.file('scripts', package='fda')) 
#> [1] ""
```

EXCEPT: These ‘scripts’ files are no longer in the package.
