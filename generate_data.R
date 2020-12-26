### data generating, assume 4 relevant dimensions
### generating y using true Z with X_j ~ N(0, 1), Z = pnorm(X)
### covariance matrix is generated as a toeplitz matrix with the
### first row being 1, 0.7, 0.5, 0.3, and 5:p being evenly
### distributed on (0.1, 0)
generate_data <- function(n, p) {
  returnList <- list()
  # generating X
  # covariance matrix
  Cov_Matrix <- toeplitz(c(1, 0.7, 0.5, 0.3, seq(0.1, 0, length.out = p-4)))
  X <- mvrnorm(n, rep(0, p), Cov_Matrix)
  # generating true Z
  Z <- apply(X, 2, pnorm)
  e <- rnorm(n) # error term
  # generating y using true Z
  y <- g1(Z[, 1]) + g2(Z[, 2]) + g3(Z[, 3]) + g4(Z[, 4]) + e
  
  returnList$y <- y
  returnList$X <- X
  returnList$e <- e # save this for oracle estimate
  returnList
}

### true functions of Z1 to Z4
g1 <- function(x) 2*exp(-2.5*x) - 2*x - 1/3
g2 <- function(x) 2*(2*x-1)^2 - 25/12
g3 <- function(x) -2*sin(2*pi*x)
g4 <- function(x) x - 1/3

g_act <- c(g1, g2, g3, g4)

### true derivatives of X1 to X4, for evaluation
f1.deriv <- function(x) -5*exp(-2.5*pnorm(x))*dnorm(x) - 2*dnorm(x)
f2.deriv <- function(x) 8*(2*pnorm(x)-1)*dnorm(x)
f3.deriv <- function(x) -4*pi*cos(2*pi*pnorm(x))*dnorm(x)
f4.deriv <- function(x) dnorm(x)

f.deriv <- c(f1.deriv, f2.deriv, f3.deriv, f4.deriv)