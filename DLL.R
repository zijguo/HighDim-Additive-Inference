### DLL.R
### Functions: implement Decorrelated Local Linear
###            Estimator for functional derivative
###            estimate in high-dimensional sparse
###            additive models as developed in the
###            paper Guo and Zhang (2019)


### cv.SAM 
### Function: Cross validated sparse additive model
###           fitting to select the best number of
###           basis functions in smoothing and the 
###           best tuning parameter lambda
### Input: X, continuous, the covariates (n by p matrix)
###        y, continuous, the outcome (n by 1 vector)
###        nfold, integer, the numbers of folds for cross validation,
###               default value is 10
###        degree, vector of integers, a sequence of number of basis 
###                functions to be considered, default is 3:4
###        lam.seq, continuous, a sequence of lambdas to be considered,
###                 usually in decreasing order to speed up
### Output: a list of
###         Z.hat, the transformed covariates in a compact support, 
###                [0, 1] in our case
###         best.sam, the sparse additive model object as in the SAM
###                  package with the best parameters
###         sigma1.sq, the consistent estimate of the variance of the
###                    noise using the mean squared errors
###         X, the input covariates
cv.SAM <- function(X, y, nfold=10, degree=3:4, lam.seq) {
  # quantile transformation of X
  Z.hat <- apply(X, 2, function(x) ecdf(x)(x))
  
  n <- nrow(Z.hat); p <- nrow(Z.hat)
  # randomly shuffle the data
  shuffle_ind <- sample(1:n)
  Z.hat <- Z.hat[shuffle_ind, ]
  y <- y[shuffle_ind]
  
  len_d <- length(degree)
  len_lam <- length(lam.seq)
  MSE <- matrix(0, len_d, 3)
  colnames(MSE) <- c("degree", "lambda", "MSE")
  
  # cross validation procedures
  # break the index into n folds
  folds <- cut(seq(1,n),breaks=nfold,labels=FALSE)
  for (i in 1:len_d) {
    mse.lam <- rep(0, len_lam)
    for (fold in 1:nfold) {
      testInd <- which(folds == fold, arr.ind = TRUE)
      # create train and test data
      test_Z <- Z.hat[testInd, ]
      test_y <- y[testInd]
      train_Z <- Z.hat[-testInd, ]
      train_y <- y[-testInd]
      
      sam.fit <- samQL(train_Z, train_y, p = degree[i], lambda = lam.seq)
      y_pred <- predict(sam.fit, newdata = test_Z)$values
      mse.lam <- mse.lam + apply(test_y-y_pred, 2, function(x) mean(x^2))/nfold
    }
    
    temp.lam <- lam.seq[which.min(mse.lam)]
    MSE[i, "degree"] <- degree[i]
    MSE[i, "lambda"] <- temp.lam
    MSE[i, "MSE"] <- min(mse.lam)
  }
  
  best.ind <- which.min(MSE[, "MSE"])
  best.degree <- MSE[best.ind, "degree"]
  best.lam <- MSE[best.ind, "lambda"]
  best.sam <- samQL(Z.hat, y, p = best.degree, lambda = best.lam)
  sigma1.sq <- mean((y-predict(best.sam, newdata = Z.hat)$values)^2)
  
  returnList <- list(Z.hat = Z.hat,
                     best.sam = best.sam,
                     sigma1.sq = sigma1.sq,
                     X = X)
  returnList
}


### predict.SAM
### Function: calculate the estimate of functions given the
###           fitted sparse additive models and the test data
### Input: sam.obj, the fitted sparse additive model in cv.SAM
###        Xt, continuous, the test covariates (nt by p matrix)
### Output: f.hat, the estimated functions (nt by p matrix)
predict.SAM <- function(sam.obj, Xt) {
  nt <- nrow(Xt)
  p <- ncol(Xt)
  X <- sam.obj$X
  best.sam <- sam.obj$best.sam
  # scale the X using quantile transformation
  emp.cdf <- apply(X, 2, function(x) ecdf(x))
  X.scale <- matrix(0, nt, p)
  for (j in 1:p) {
    X.scale[, j] <- emp.cdf[[j]](Xt[, j])
  }
  
  # calculate each estimated function
  f.hat <- matrix(0, nt, p)
  for (j in 1:p) {
    bspline <- ns(X.scale[, j], df = best.sam$p, knots = best.sam$knots[, j],
                  Boundary.knots = best.sam$Boundary.knots[, j]) # B-spline matrix for prediction
    ind.j <- (j-1)*best.sam$p + c(1:best.sam$p) # index of the target component function in each loop
    f.hat[, j] <- bspline %*% best.sam$w[ind.j]
  }
  return(f.hat)
}


### LASSO
### Function: computes the LASSO estimator
### -If lambda is given, use glmnet and standard Lasso
### -If lambda is set to the character string "CV", then glmnet with 
### lambda selected by cross-validation is used
### -If lambda is not given or is set to NULL, use square root Lasso
Lasso <- function(X, y, lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)
  
  htheta <- if (is.null(lambda)) {
    lambda <- sqrt(qnorm(1 - (0.1 / p)) / n)
    outLas <- slim(X, y, lambda = lambda, method = "lq", q = 2,
                   verbose = FALSE)
    # Objective : sqrt(RSS/n) + lambda * penalty
    c(as.vector(outLas$intercept), as.vector(outLas$beta))
  } else if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "CV.min") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  } else if (lambda == "scalreg") {
    Xc <- if (intercept) {
      cbind(rep(1, n), X)
    } else {
      X
    }
    outLas <- scalreg(Xc, y)
    # return object
    if (intercept) {
      outLas$coefficients
    } else {
      # add a coefficient for the (not estimated) intercept b/c of implementation
      c(0, outLas$coefficients)
    }
  } else {
    outLas <- glmnet(X, y, family = "gaussian", alpha = 1,
                     intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = lambda))
  }
  
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}


### Initialization.step
### Function: computes the initial LASSO estimator and quantities 
### based thereon
### lambda is set as "CV"
Initialization.step <- function(X, y, lambda = NULL, intercept = FALSE) {
  n <- nrow(X)
  # col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
  col.norm <- 1 / sqrt((1 / n) * diag(t(X)%*%X))
  Xnor <- X %*% diag(col.norm)
  
  ### Call Lasso
  htheta <- Lasso(Xnor, y, lambda = lambda, intercept = intercept)
  
  ### Calculate return quantities
  if (intercept == TRUE) {
    Xb <- cbind(rep(1, n), Xnor)
    col.norm <- c(1, col.norm)
  } else {
    Xb <- Xnor
  }
  sparsity <- sum(abs(htheta) > 0.001)
  sd.est <- sqrt(sum((y - Xb %*% htheta)^2) / n)
  htheta <- htheta * col.norm
  returnList <- list("lasso.est" = htheta,
                     "sigma" = sd.est,
                     "sparsity" = sparsity)
  return(returnList)
}


### DLL
### Function: implement the Decorrelated Local Linear 
###           estimator based on the initial estimate
###           of functions using cv.SAM
### input: X, continuous, the covariates (n by p matrix)
###        y, continuous, the outcome (n by 1 vector)
###        sam.obj, the sparse additive model object list
###                 obtained from cv.SAM
###        x.eval, vector, the evaluation points to estimate the derivatives
###        ind, column index of the function of interest in the covariates
###        h, optional, a single value or a vector of the same length as 
###           x.eval, the bandwidth for DLL estimator
###        dat.split, logical, implement the data swapping method in the
###                   paper, the default value is TRUE
###                   
### output: a list of
###         x.eval, vector, the evaluation points to estimate the derivatives
###         est, the point estimates of derivatives at x.eval
###         est.se, the standard errors of the DLL estimates at x.eval
###         Sn, the Sn estimated at x.eval
###         l.na.index, a vector of logical values indicating where NA appears
###         for estimating l due to the 0 denominator
DLL <- function(X, y, sam.obj, x.eval, ind, h=NULL) {
  n <- nrow(X)
  p <- ncol(X)
  n.eval <- length(x.eval)
  
  f.hat <- predict.SAM(sam.obj, X)
  sigma1.sq <- sam.obj$sigma1.sq # get sigma1.sq for variance estimation
  
  # calculate R by detracting nuisance function
  R.hat <- y - apply(f.hat[,-ind], 1, sum)
  
  if (is.null(h)) {
    # h <- npregbw(R.hat~X[, ind], regtype="ll")$bw
    # h <- n^(-0.2)
    h <- lpbwselect(R.hat, X[, ind], eval = x.eval, p = 1, bwselect = "ce-dpi", deriv = 1,
                    kernel = "uni")$bws[, "h"]
    # h <- regCVBwSelC(X[, ind], R.hat, deg = 1, kernel = SqK)
    # h <- thumbBw(X[, ind], R.hat, deg = 1, kernel = SqK)
  }
  
  # check for single bandwidth or a vector
  if (length(h)==1) h <- rep(h, n.eval)
  
  # LASSO estimator to get delta.hat and mu.hat
  Ini.Est <- Initialization.step(X[,-ind], X[,ind], lambda = "CV.min", intercept = TRUE)
  
  est <- est.se <- Sn <- rep(0, n.eval) 
  for (j in 1:n.eval) {
    delta <- matrix(X[, ind]) - cbind(1, matrix(X[, -ind], n, p-1)) %*% Ini.Est$lasso.est
    mu <- rep(x.eval[j], n) - cbind(1, matrix(X[, -ind], n, p-1)) %*% Ini.Est$lasso.est
    
    ### calculate uniform Kernels
    Kh <- ifelse(abs(X[, ind]-rep(x.eval[j], n)) <= h[j], yes = 1/h[j], no = 0)
    
    ### estimate Di1
    l <- rep(0, n)
    for (i in 1:n) {
      wt <- sum(ifelse(abs(delta - rep(mu[i], n)) <= h[j], yes = 1, no = 0))
      l[i] <- sum((delta - rep(mu[i], n))*ifelse(abs(delta - rep(mu[i], n)) <= h[j], yes = 1, no = 0))/wt
    }
    # if abs(delta-mu) is not in the bandwidth, set l to be zero
    l.na.ind <- is.na(l)
    l[l.na.ind] <- 0
    
    D.tilde <- (X[, ind] - rep(x.eval[j], n)) - l
    D.hat <- D.tilde - sum(Kh*D.tilde)/sum(Kh)
    ### calculate Sn
    Sn.j <- 1/n*sum(D.hat*(X[, ind]-rep(x.eval[j], n))*Kh)
    
    est[j] <- 1/(n*Sn.j)*sum(D.hat*R.hat*Kh)
    Sn[j] <- Sn.j
    est.se[j] <- sqrt(sigma1.sq/(n^2*Sn.j^2)*sum(D.hat^2*Kh^2))
  }
  
  returnList <- list(x.eval = x.eval,
                     est = est,
                     est.se = est.se,
                     Sn = Sn,
                     l.na.ind = l.na.ind)
  returnList
}
