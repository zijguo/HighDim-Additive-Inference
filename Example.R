### required packages
library(SAM) # sparse additive model
library(splines) # generating spline matrix
library(glmnet) # lasso estimator
library(MASS) # generating multivariate standard normal distribution
library(nprobust) # functions of bandwidth selection

n <- 500
p <- 750

# generate data
x.eval <- c(-1.75, 0.1, 0.3, 0.5, 0.75)
set.seed(1234)
data_list <- generate_data(n, p)
y <- data_list$y
X <- data_list$X
e <- data_list$e

# sparse additive model
out.sam <- cv.SAM(X, y, degree = 3, lam.seq = 0.5^seq(0, 10, length.out = 200))

# the true derivatives
f1.deriv(x.eval)

# DLL
head(DLL(X, y, out.sam, x.eval, ind = 1), n=4)

# plug-in
loclin.plug(X, y, out.sam, x.eval, ind = 1)

# oracle
loclin.orac(X, y, x.eval, g_act, pnorm, e, ind = 1)
