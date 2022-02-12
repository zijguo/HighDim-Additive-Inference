
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DLL

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/DLL)](https://CRAN.R-project.org/package=DLL)
<!-- badges: end -->

The goal of DLL is to implement the Decorrelated Local Linear estimator
proposed in &lt;arxiv:1907.12732&gt;. It constructs the confidence
interval for the derivative of the function of interest under the
high-dimensional sparse additive model.

## Installation

You can install the released version of DLL from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("DLL")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(DLL)
library(MASS)
# evaluation points
d0 = c(-0.5,0.25)

f = function(x) 1.5*sin(x)
f.deriv = function(x) 1.5*cos(x)
g1 = function(x) 2*exp(-x/2)
g2 = function(x) (x-1)^2 - 25/12
g3 = function(x) x - 1/3
g4 = function(x) 0.75*x
g5 = function(x) 0.5*x


# sample size and dimension of X
n = 500
p = 500

# covariance structure of D and X
Cov_Matrix = toeplitz(c(1, 0.7, 0.5, 0.3, seq(0.1, 0, length.out = p-3)))

set.seed(123)
# X represents the (D,X) here
X = mvrnorm(n,rep(-0.25,p+1),Sigma = Cov_Matrix)
e = rnorm(n,sd=1)
# generating response
y = f(X[,1]) + g1(X[,2]) + g2(X[,3]) + g3(X[,4]) + g4(X[,5]) + g5(X[,6]) + e

### DLL inference
DLL.model = DLL(X=X, y=y, D.ind = 1, d0 = d0)
```

true values

``` r
f.deriv(d0)
#> [1] 1.316374 1.453369
```

point estimates

``` r
DLL.model$est
#>            f1
#> -0.5 1.258581
#> 0.25 1.659544
```

standard errors

``` r
DLL.model$est.se
#>             f1
#> -0.5 0.3911074
#> 0.25 0.4301377
```

confidence interval

``` r
DLL.model$CI
#> $lower
#>             f1
#> -0.5 0.4920249
#> 0.25 0.8164900
#> 
#> $upper
#>            f1
#> -0.5 2.025138
#> 0.25 2.502599
```
