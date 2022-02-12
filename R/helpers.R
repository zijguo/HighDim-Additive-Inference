### cv.SAM
### Function: Cross-validated Sparse Additive Model
###           fitting to select the best number of
###           basis functions in smoothing and the
###           best tuning parameter lambda
### Input: X, continuous, the covariates (D,X)
###        y, continuous, the outcome
###        kfold, integer, the numbers of folds for cross validation,
###               default value is 5
###        degree, vector of integers, a sequence of number of basis
###                functions to be considered, default is 3:6
###        lam.seq, continuous, a sequence of lambdas to be considered,
###                 usually in decreasing order to speed up, default to
###                 exp(seq(log(1),log(5e-3),length = 100)) if NULL
###        quant.trans, to conduct quantile transformation or not
### Output: sam.final, the sparse additive model object as in the SAM
###                  package with the best parameters
###         sigma1.sq, the consistent estimate of the variance of the
###                    noise using the mean squared errors
###         X, input covariates (D,X)
###         Z.hat, the transformed covariates (D,X) in a compact support,
###                [0, 1] if using quantile transformation in our case, same to
###                X if quant.trans=FALSE
cv.SAM = function(X, y, kfold=5, degree=3:6, lam.seq, quant.trans) {
  n = nrow(X); p = ncol(X)
  # quantile transformation of X
  if (quant.trans) {
    Z.hat = apply(X, 2, function(x) ecdf(x)(x))
  } else {
    Z.hat = X
  }

  # default lambda sequence, approximately from 0.001 to 1
  if (is.null(lam.seq)) lam.seq = exp(seq(log(1),log(5e-3),length = 100))

  # cross validation
  len_d = length(degree)
  len_lam = length(lam.seq)
  MSE = matrix(0, len_d, 3)
  colnames(MSE) = c("degree", "lambda", "MSE")
  # break the index into k folds
  folds = cut(seq(1,n),breaks=kfold,labels=FALSE)
  for (i in 1:len_d) {
    mse.lam = rep(0, len_lam)
    for (fold in 1:kfold) {
      testInd = which(folds == fold, arr.ind = TRUE)
      # create train and test data
      Z.test = Z.hat[testInd, ]
      y.test = y[testInd]
      Z.train = Z.hat[-testInd, ]
      y.train = y[-testInd]
      # provide a sequence of lambda to speed up
      sam.fit = samQL(Z.train, y.train, p = degree[i], lambda = lam.seq)
      y.pred = predict(sam.fit, newdata = Z.test)$values
      mse.lam = mse.lam + apply(y.test-y.pred, 2, function(x) mean(x^2))/kfold
    }

    temp.lam = lam.seq[which.min(mse.lam)]
    MSE[i, "degree"] = degree[i]
    MSE[i, "lambda"] = temp.lam
    MSE[i, "MSE"] = min(mse.lam)
  }

  best.ind = which.min(MSE[, "MSE"])
  best.degree = MSE[best.ind, "degree"]
  best.lam = MSE[best.ind, "lambda"]
  sam.final = samQL(Z.hat, y, p = best.degree, lambda = best.lam)
  sigma1.sq = mean((y-predict(sam.final, newdata = Z.hat)$values)^2)

  returnList = list(sam.final = sam.final,
                    sigma1.sq = sigma1.sq,
                    X = X,
                    Z.hat = Z.hat,
                    quant.trans = quant.trans)
  returnList
}


### predict.SAM
### Function: Calculate the estimate of functions given the
###           fitted sparse additive models and the test data
### Input: sam.obj, fitted sparse additive model from cv.SAM
###        Xt, continuous, the test covariates (nt by p matrix)
### Output: f.hat, estimated functions including an intercept
###                (nt by p+1 matrix)
predict.SAM = function(sam.obj, Xt) {
  quant.trans = sam.obj$quant.trans
  nt = nrow(Xt)
  p = ncol(Xt)
  X = sam.obj$X
  sam.final = sam.obj$sam.final
  # scale the X using quantile transformation or not
  if (quant.trans) {
    emp.cdf = apply(X, 2, function(x) ecdf(x))
    X.qt = matrix(0, nt, p)
    for (j in 1:p) {
      X.qt[, j] = emp.cdf[[j]](Xt[, j])
    }
  } else {
    X.qt = Xt
  }
  # max-min transformation default in the SAM package
  X.min.rep = matrix(rep(sam.final$X.min,nt),nrow=nt,byrow=T)
  X.ran.rep = matrix(rep(sam.final$X.ran,nt),nrow=nt,byrow=T)

  X.scale = (X.qt-X.min.rep)/X.ran.rep
  X.scale = pmax(X.scale,0)
  X.scale = pmin(X.scale,1)


  # calculate each estimated function
  f.hat = matrix(0, nt, p)
  for (j in 1:p) {
    bspline = ns(X.scale[, j], df = sam.final$p, knots = sam.final$knots[, j],
                 Boundary.knots = sam.final$Boundary.knots[, j]) # B-spline matrix for prediction
    ind.j = (j-1)*sam.final$p + c(1:sam.final$p) # index of the target component function coefficients
    f.hat[, j] = bspline %*% sam.final$w[ind.j]
  }
  # add the intercept
  f.hat = cbind(rep(sam.final$intercept, nt), f.hat)
  colnames(f.hat) = c("Intercept", paste("X", 1:p, sep = ""))

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

  htheta <- if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "CV.min") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  }
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}
# Lasso = function(X, y, lambda = NULL, intercept = TRUE) {
#   p = ncol(X)
#   n = nrow(X)
#
#   htheta = if (is.null(lambda)) {
#     lambda = sqrt(qnorm(1 - (0.1 / p)) / n)
#     outLas = slim(X, y, lambda = lambda, method = "lq", q = 2,
#                   verbose = FALSE)
#     # Objective : sqrt(RSS/n) + lambda * penalty
#     c(as.vector(outLas$intercept), as.vector(outLas$beta))
#   } else if (lambda == "CV") {
#     outLas = cv.glmnet(X, y, family = "gaussian", alpha = 1,
#                        intercept = intercept)
#     # Objective : 1/2 * RSS/n + lambda * penalty
#     as.vector(coef(outLas, s = outLas$lambda.1se))
#   } else if (lambda == "CV.min") {
#     outLas = cv.glmnet(X, y, family = "gaussian", alpha = 1,
#                        intercept = intercept)
#     # Objective : 1/2 * RSS/n + lambda * penalty
#     as.vector(coef(outLas, s = outLas$lambda.min))
#   } else if (lambda == "scalreg") {
#     Xc = if (intercept) {
#       cbind(rep(1, n), X)
#     } else {
#       X
#     }
#     outLas = scalreg(Xc, y)
#     # return object
#     if (intercept) {
#       outLas$coefficients
#     } else {
#       # add a coefficient for the (not estimated) intercept b/c of implementation
#       c(0, outLas$coefficients)
#     }
#   } else {
#     outLas = glmnet(X, y, family = "gaussian", alpha = 1,
#                     intercept = intercept)
#     # Objective : 1/2 * RSS/n + lambda * penalty
#     as.vector(coef(outLas, s = lambda))
#   }
#
#   if (intercept == TRUE) {
#     return(htheta)
#   } else {
#     return(htheta[2:(p+1)])
#   }
# }



### Initialization.step
### Function: Computes the initial LASSO estimator and quantities based thereon
### lambda is set as "CV"
Initialization.step = function(X, y, lambda = NULL, intercept = FALSE) {
  n = nrow(X)
  # col.norm = 1 / sqrt((1 / n) * diag(t(X) %*% X))
  col.norm = 1 / sqrt((1 / n) * diag(t(X)%*%X))
  Xnor = X %*% diag(col.norm)

  ### Call Lasso
  htheta = Lasso(Xnor, y, lambda = lambda, intercept = intercept)

  ### Calculate return quantities
  if (intercept == TRUE) {
    Xb = cbind(rep(1, n), Xnor)
    col.norm = c(1, col.norm)
  } else {
    Xb = Xnor
  }
  sparsity = sum(abs(htheta) > 0.001)
  sd.est = sqrt(sum((y - Xb %*% htheta)^2) / n)
  htheta = htheta * col.norm
  returnList = list("lasso.est" = htheta,
                    "sigma" = sd.est,
                    "sparsity" = sparsity)
  return(returnList)
}
