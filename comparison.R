### comparison.R
### Functions: use the plug-in estimator and the
###            oracle estimator for comparison

### loclin.plug
### Function: the plug-in local linear estimator,
###           subtract y by f_2 to f_p estimated 
###           from sparse additive model to get
###           R.hat, fit on R.hat~X_1
### The input and output are similar to DLL
loclin.plug <- function(X, y, sam.obj, x.eval, ind, h=NULL) {
  n <- nrow(X)
  p <- ncol(X)
  n.eval <- length(x.eval)
  
  f.hat <- predict.SAM(sam.obj, X)
  # calculate R by detracting nuisance function
  R.hat <- y - apply(f.hat[,-ind], 1, sum)
  
  ## calculate sigma1.sq, get estimates of local linear estimator at each data point
  h.sigma1 <- lpbwselect(R.hat, X[, ind], eval = X[, ind], p = 1, bwselect = "mse-dpi",
                         kernel = "uni", bwregul = 0)$bws[, "h"]
  # h.sigma1 <- n^(-0.2)
  # check for single bandwidth or a vector
  if (length(h.sigma1)==1) h.sigma1 <- rep(h.sigma1, n)
  # estimate of sigma1.sq for local linear estimator
  res <- rep(0, n) # residuals
  for (j in 1:n) {
    ### calculate uniform Kernels
    Kh <- ifelse(abs(X[, ind]-rep(X[j, ind], n)) <= h.sigma1[j], yes = 1/h.sigma1[j], no = 0)
    D.tilde <- (X[, ind] - rep(X[j, ind], n))
    # calculate sigma1.square
    W <- diag(Kh)
    XX <- cbind(1, D.tilde)
    est <- t(matrix(c(1,0)))%*%solve(t(XX)%*%W%*%XX)%*%t(XX)%*%W%*%R.hat
    res[j] <- R.hat[j] - est
  }
  sigma1.sq <- mean(res^2) # estimate of sigma1.sq
  print(sigma1.sq)
  
  if (is.null(h)) {
    # h <- npregbw(R.hat~X[, ind], regtype="ll")$bw
    # h <- n^(-0.2)
    # h <- lpbwselect(R.hat, X[, ind], eval = x.eval, p = 1, bwselect = "mse-dpi",
    #                 kernel = "uni")$bws[, "h"]
    # h <- regCVBwSelC(X[, ind], R.hat, deg = 1, kernel = SqK)
    h <- thumbBw(X[, ind], R.hat, deg = 1, kernel = SqK)
  }
  # check for single bandwidth or a vector
  if (length(h)==1) h <- rep(h, n.eval)
  
  est.plug <- est.plug.se <- Sn.plug <- rep(0, n.eval)
  for (j in 1:n.eval) {
    ### calculate uniform Kernels
    Kh <- ifelse(abs(X[, ind]-rep(x.eval[j], n)) <= h[j], yes = 1/h[j], no = 0)
    D.tilde <- (X[, ind] - rep(x.eval[j], n))
    D.hat <- D.tilde - sum(Kh*D.tilde)/sum(Kh)
    ### calculate Sn
    Sn.j <- 1/n*sum(D.hat*(X[, ind]-rep(x.eval[j], n))*Kh)
    
    est.plug[j] <- 1/(n*Sn.j)*sum(D.hat*R.hat*Kh)
    est.plug.se[j] <- sqrt(sigma1.sq/(n^2*Sn.j^2)*sum(D.hat^2*Kh^2))
    Sn.plug[j] <- Sn.j
  }
  
  returnList <- list(x.eval = x.eval,
                     est.plug = est.plug,
                     est.plug.se = est.plug.se,
                     Sn.plug = Sn.plug)
  returnList
}


### loclin.orac
### Function: the oracle local linear estimator,
###           subtract y by true f_2 to f_p to 
###           get R.true, fit on R.true~X_1
### The input and output are similar to DLL
### g_act, a list of true functions at each dimension
###        used to construct R.true
### true.cdf, the true CDF of X_1
### e, continuous, the error term (n by 1 vector), used
###    to construct R.true
loclin.orac <- function(X, y, sam.obj, x.eval, g_act, true.cdf, e, ind, h=NULL) {
  n <- nrow(X)
  p <- ncol(X)
  n.eval <- length(x.eval)
  
  # calculate R by detracting nuisance function
  f.hat <- predict.SAM(sam.obj, X)
  R.hat <- y - apply(f.hat[,-ind], 1, sum)
  # calculate R.true by detracting true nuisance function
  R.true <- g_act[[ind]](true.cdf(X[, ind])) + e
  
  ## calculate sigma1.sq, get estimates of local linear estimator at each data point
  h.sigma1 <- lpbwselect(R.true, X[, ind], eval = X[, ind], p = 1, bwselect = "mse-dpi",
                         kernel = "uni", bwregul = 0)$bws[, "h"]
  # h.sigma1 <- n^(-0.2)
  # check for single bandwidth or a vector
  if (length(h.sigma1)==1) h.sigma1 <- rep(h.sigma1, n)
  # estimate of sigma1.sq for local linear estimator
  res <- rep(0, n) # residuals
  for (j in 1:n) {
    ### calculate uniform Kernels
    Kh <- ifelse(abs(X[, ind]-rep(X[j, ind], n)) <= h.sigma1[j], yes = 1/h.sigma1[j], no = 0)
    D.tilde <- (X[, ind] - rep(X[j, ind], n))
    # calculate sigma1.square
    W <- diag(Kh)
    XX <- cbind(1, D.tilde)
    est <- t(matrix(c(1,0)))%*%solve(t(XX)%*%W%*%XX)%*%t(XX)%*%W%*%R.true
    res[j] <- R.true[j] - est
  }
  sigma1.sq <- mean(res^2) # estimate of sigma1.sq
  print(sigma1.sq)
  
  # use the same bandwidth for all estimators
  if (is.null(h)) {
    # h <- npregbw(R.true~X[, ind], regtype="ll")$bw
    # h <- n^(-0.2)
    # h <- lpbwselect(R.true, X[, ind], eval = x.eval, p = 1, bwselect = "mse-dpi",
    #                 kernel = "uni")$bws[, "h"]
    # h <- regCVBwSelC(X[, ind], R.true, deg = 1, kernel = SqK)
    h <- thumbBw(X[, ind], R.hat, deg = 1, kernel = SqK) # use the same bandwidth
  }
  # check for single bandwidth or a vector
  if (length(h)==1) h <- rep(h, n.eval)
  
  est.orac <- est.orac.se <- Sn.orac <- rep(0, n.eval)
  for (j in 1:n.eval) {
    ### calculate uniform Kernels
    Kh <- ifelse(abs(X[, ind]-rep(x.eval[j], n)) <= h[j], yes = 1/h[j], no = 0)
    D.tilde <- (X[, ind] - rep(x.eval[j], n))
    D.hat <- D.tilde - sum(Kh*D.tilde)/sum(Kh)
    ### calculate Sn
    Sn.j <- 1/n*sum(D.hat*(X[, ind]-rep(x.eval[j], n))*Kh)
    
    est.orac[j] <- 1/(n*Sn.j)*sum(D.hat*R.true*Kh)
    est.orac.se[j] <- sqrt(sigma1.sq/(n^2*Sn.j^2)*sum(D.hat^2*Kh^2))
    Sn.orac[j] <- Sn.j
  }
  
  returnList <- list(x.eval = x.eval,
                     est.orac = est.orac,
                     est.orac.se = est.orac.se,
                     Sn.orac = Sn.orac)
  returnList
}
