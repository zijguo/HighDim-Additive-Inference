### source DLL.R and comparison.R before running the code


### required packages
library(SAM) # sparse additive models
library(splines) # generating spline matrix
library(glmnet) # lasso estimator
library(MASS) # generating multivariate standard normal distribution
library(nprobust) # estimator of bandwidth for sigma1.sq
library(locpol) # estimate of bandwidth for DLL
library(np) # bandwidth selection for DLL

mean.x <- 0.25; sd.x <- 1
p.x <- function(x) pnorm(x, mean = mean.x, sd = sd.x) # CDF of X
d.x <- function(x) dnorm(x, mean = mean.x, sd = sd.x) # PDF of X

generate_data <- function(n, p, g_act) {
  returnList <- list()
  # generating X
  # covariance matrix
  Cov_Matrix <- toeplitz(c(1, 0.7, 0.5, 0.3, seq(0.1, 0, length.out = p-4)))
  X <- mvrnorm(n, rep(mean.x, p), Cov_Matrix)
  # generating true Z
  Z <- apply(X, 2, p.x)
  e <- rnorm(n) # error term
  y <- e
  # generating y using true Z
  # approximately sparse and number of g should not being tested jointly
  if(approx_sparse==1) {
    if (one_non_linear==0) {
      y <- y + g1(Z[,1]) + g2(Z[,2]) + g3(Z[,3]) + g4(Z[,4])
      for (i in 5:p) {
        y <- y + (i-1)^(-1)*Z[, i]
      }
    } else if (one_non_linear==1) {
      y <- y + g1(Z[,1])
      for (i in 2:p) {
        y <- y + (i-1)^(-1)*Z[, i]
      }
    }
  } else {
    for (i in 1:length(g_act)) {
      y <- y + g_act[[i]](Z[, i]) # general case
    }
  }
  returnList$y <- y
  returnList$X <- X
  returnList$e <- e # save this for oracle estimate
  returnList
}

### true functions of Z1 to Z4
g1<- function(z) 2*exp(-2.5*z) - 2*z - 1/3
g2 <- function(z) 2*(2*z-1)^2 - 25/12
g3 <- function(z) -sin(2*pi*z)
g4 <- function(z) z - 1/3

## this is adding some linear and non-linear terms
g5 <- function(z) 0.5*z
g6 <- function(z) 0.4*z
g7 <- function(z) 0.3*z
g8 <- function(z) 0.2*z
g9 <- function(z) 0.1*z
g10 <- function(z) 0.1*sin(2*pi*z)
g11 <- function(z) 0.2*cos(2*pi*z)
g12 <- function(z) 0.3*(sin(2*pi*z))^2
g13 <- function(z) 0.4*(cos(2*pi*z))^3
g14 <- function(z) 0.5*(sin(2*pi*z))^3
g15 <- function(z) z/(1+z)

### true derivatives of X1 to X4, for evaluation
f1.deriv <- function(x) -5*exp(-2.5*p.x(x))*d.x(x) - 2*d.x(x)
f2.deriv <- function(x) 8*(2*p.x(x)-1)*d.x(x)
f3.deriv <- function(x) -2*pi*cos(2*pi*p.x(x))*d.x(x)
f4.deriv <- function(x) d.x(x)

f.deriv <- c(f1.deriv, f2.deriv, f3.deriv, f4.deriv)

approx_sparse <- 0 # approximately sparse if 1
one_non_linear <- 0 # only the first dimension is non-linear, use with approx_sparse=1
gnum_setting <- 0 # use 15 activation functions if 1, not with approx_sparse

# gnum_setting=1 means more sparsity(15 functions)
if (gnum_setting==1) {
  g_act <- c(g1, g2, g3, g4, g5, g6, g7, g8,
             g9, g10, g11, g12, g13, g14, g15)
} else if (gnum_setting==0) {
  g_act <- c(g1, g2, g3, g4)
}



### Example
# generate data
n <- 500 # sample size
p <- 750 # number of dimensions

# evaluation points
x.eval <- c(-1, -0.5, 0, 0.5, 1)

set.seed(2021)
data_list <- generate_data(n, p, g_act)
y <- data_list$y
X <- data_list$X
e <- data_list$e

# initial estimator, the sparse additive model fitting
out.sam <- cv.SAM(X, y, degree=3)
out.sam$sigma1.sq # the estimate of noise level

# point estimator of f.hat
f.hat <- predict.SAM(out.sam, X)

# true value of f^{'}_1 at x.eval
f1.deriv(x.eval)

# point estimate of DLL at x.eval of the first dimension
DLL.out <- DLL(X, y, out.sam, x.eval, ind=1)
DLL.out$est # point estimate
DLL.out$est.se # the standard error

# point estimate of Plug-in at x.eval of the first dimension
plug.out <- loclin.plug(X, y, out.sam, x.eval, ind=1)
plug.out$est.plug # point estimate
plug.out$est.plug.se # the standard error

# point estimate of Oracle at x.eval of the first dimension
orac.out <- loclin.orac(X, y, out.sam, x.eval, g_act, p.x, e, ind=1)
orac.out$est.orac # point estimate
orac.out$est.orac.se # the standard error
