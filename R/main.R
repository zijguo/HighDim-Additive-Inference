#' DLL
#' @description It constructs the Decorrelated Local Linear estimator and estimates its standard error.
#' It further constructs the confidence interval for the derivative of the function of interest.
#' @param X the covariates matrix, of dimension \eqn{n \times p}
#' @param y the outcome vector, of length \eqn{n}
#' @param D.ind the column index(es) of X, indicating the index(es) of the variable(s) of interest.
#' It can be a scalar or a vector. If vector, then do inference for each index of the sequence.
#' @param d0 evaluation points for derivative estimation. It can be scalar or vector.
#' @param h bandwidth, computed by Rule of Thumb from the package ``locpol'' if not provided.
#' @param lam.seq a sequence of tuning parameters considered in fitting the sparse additive model. Cross validation is used to choose the best one.
#' If not provided(default), the default sequence ranges from 5e-3 to 1 with the length of 100. If provided, the sequence needs to be in a decreasing order
#' for the reason of computation efficiency.
#' @param treatment.SAM Whether a sparse additive model is used for fitting the treatment model? If `False`(default), Lasso with cross validation is used to fit the treatment model.
#' Default is `FALSE`
#' @param data.swap Whether data swapping is conducted or not? Default is `FALSE`
#' @param quant.trans Whether quantile transformation is conducted or not? Default is `FALSE`
#' @param alpha the significance level. Default is 0.05
#'
#' @return
#' \item{est}{point estimates of the function derivative}
#' \item{est.se}{estimated standard errors of est}
#' \item{CI}{list of lower and upper bounds of confidence intervals}
#' \item{d0}{evaluation points}
#' \item{bw.save}{selected bandwidth at each element of d0}
#' \item{sigma1.sq}{estimated variance of the error term in the outcome model}
#' @export
#' @import SAM glmnet locpol
#' @importFrom stats coef ecdf predict toeplitz qnorm
#' @importFrom splines ns
#' @examples
#' # evaluation points
#' d0 = c(-0.5,0.25)
#' f = function(x) 1.5*sin(x)
#' f.deriv = function(x) 1.5*cos(x)
#' g1 = function(x) 2*exp(-x/2)
#' g2 = function(x) (x-1)^2 - 25/12
#' g3 = function(x) x - 1/3
#' g4 = function(x) 0.75*x
#' g5 = function(x) 0.5*x
#' # sample size and dimension of X
#' n = 200
#' p = 100
#' # covariance structure of D and X
#' Cov_Matrix = toeplitz(c(1, 0.7, 0.5, 0.3, seq(0.1, 0, length.out = p-3)))
#' set.seed(123)
#' # X represents the (D,X) here
#' X = MASS::mvrnorm(n,rep(-0.25,p+1),Sigma = Cov_Matrix)
#' e = rnorm(n,sd=1)
#' # generating response
#' y = f(X[,1]) + g1(X[,2]) + g2(X[,3]) + g3(X[,4]) + g4(X[,5]) + g5(X[,6]) + e
#' ### DLL inference
#' DLL.model = DLL(X=X, y=y, D.ind = 1, d0 = d0)
#' # true values
#' f.deriv(d0)
#' # point estimates
#' DLL.model$est
#' # standard errors
#' DLL.model$est.se
#' # confidence interval
#' DLL.model$CI
DLL = function(X, y, D.ind, d0, h=NULL, lam.seq = NULL, treatment.SAM=FALSE, data.swap=FALSE, quant.trans=FALSE,alpha=0.05) {
  n = nrow(X)
  p = ncol(X)
  # center X and d0 before fitting
  d0.c = d0 - mean(X[,D.ind])
  X = X - matrix(apply(X,2,mean),n,p,byrow = TRUE)


  n.d0 = length(d0)
  if (data.swap) {
    a.ind = 1:round(n/2); b.ind = setdiff(1:n,a.ind)
    n.a = length(a.ind); n.b = length(b.ind)
  }

  # fit sparse additive model
  if (data.swap) {
    sam.a = cv.SAM(X[a.ind,],y[a.ind],lam.seq = lam.seq,quant.trans=quant.trans)
    sam.b = cv.SAM(X[b.ind,],y[b.ind],lam.seq = lam.seq,quant.trans=quant.trans)
    f.hat.a = predict.SAM(sam.b,X[a.ind,])
    f.hat.b = predict.SAM(sam.a,X[b.ind,])
    f.hat = rbind(f.hat.a,f.hat.b)
  } else {
    sam.model = cv.SAM(X,y,lam.seq = lam.seq,quant.trans=quant.trans)
    f.hat = predict.SAM(sam.model,X)
  }

  sigma1.sq = mean((y-apply(f.hat,1,sum))^2) # get sigma1.sq for variance estimation

  # calculate point estimator and standard error
  est = est.se = matrix(NA,length(d0),length(D.ind))
  colnames(est) = colnames(est.se) = paste("f",as.character(D.ind),sep="")
  rownames(est) = rownames(est.se) = d0
  bw.save = rep(list(NA),length(D.ind))

  for (d.ind in 1:length(D.ind)) {
    # calculate R by detracting nuisance function
    R.hat = y - apply(f.hat[,-(D.ind[d.ind]+1)], 1, sum) # D.ind[d.ind] plus 1 to skip the intercept
    h = suppressWarnings(thumbBw(X[,D.ind[d.ind]],R.hat,deg=1,kernel=SqK))
    # h = lpbwselect(R.hat,X[,D.ind[d.ind]],eval=d0.c,deriv=1,kernel="uni",bwselect="mse-dpi")$bws[,"h"]
    # h = lpbwselect(R.hat,X[,D.ind[d.ind]],eval=d0.c,deriv=1,kernel="uni",bwselect="ce-rot")$bws[,"h"]
    # check for single bandwidth or a vector
    if (length(h)==1) h = rep(h, n.d0)
    bw.save[[d.ind]] = h

    if (data.swap) {
      # LASSO estimator to get delta.hat and mu.hat
      if (!treatment.SAM) {
        gamma.a = Initialization.step(X[a.ind,-D.ind[d.ind]], X[a.ind,D.ind[d.ind]], lambda = "CV.min", intercept = TRUE)$lasso.est
        gamma.b = Initialization.step(X[b.ind,-D.ind[d.ind]], X[b.ind,D.ind[d.ind]], lambda = "CV.min", intercept = TRUE)$lasso.est
      } else {
        # sam model for treatment model D ~ X
        sam.D.a = cv.SAM(X[a.ind,-D.ind[d.ind]],X[a.ind,D.ind[d.ind]],lam.seq = lam.seq,quant.trans = FALSE)
        sam.D.b = cv.SAM(X[b.ind,-D.ind[d.ind]],X[b.ind,D.ind[d.ind]],lam.seq = lam.seq,quant.trans = FALSE)
        D.hat.a = predict(sam.D.b$sam.final,newdata=X[a.ind,-D.ind[d.ind]])$values
        D.hat.b = predict(sam.D.a$sam.final,newdata=X[b.ind,-D.ind[d.ind]])$values
      }

      for (j in 1:n.d0) {
        # calculate mu and delta
        if (!treatment.SAM) {
          mu.a = d0.c[j] - cbind(1,matrix(X[a.ind,-D.ind[d.ind]],n.a,p-1))%*%gamma.b
          mu.b = d0.c[j] - cbind(1,matrix(X[b.ind,-D.ind[d.ind]],n.b,p-1))%*%gamma.a
          delta.a = matrix(X[a.ind,D.ind[d.ind]],n.a,1) - cbind(1,matrix(X[a.ind,-D.ind[d.ind]],n.a,p-1))%*%gamma.b
          delta.b = matrix(X[b.ind,D.ind[d.ind]],n.b,1) - cbind(1,matrix(X[b.ind,-D.ind[d.ind]],n.b,p-1))%*%gamma.a
          delta = c(delta.a,delta.b)
        } else {
          mu.a = d0.c[j] - D.hat.a
          mu.b = d0.c[j] - D.hat.b
          delta.a = matrix(X[a.ind,D.ind[d.ind]],n.a,1) - D.hat.a
          delta.b = matrix(X[b.ind,D.ind[d.ind]],n.b,1) - D.hat.b
          delta = c(delta.a,delta.b)
        }
        # calculate l
        l.a = rep(0,n.a); l.b = rep(0,n.b)
        for (i in 1:n.a) {
          # weight
          w.a = ifelse(abs(delta.a-mu.a[i])<=h[j],yes=1,no=0)
          l.a[i] = sum((delta.a-mu.a[i])*w.a)/sum(w.a)
          l.a[is.na(l.a)] = 0 # remove NAs
        }
        for (i in 1:n.b) {
          # weight
          w.b = ifelse(abs(delta.b-mu.b[i])<=h[j],yes=1,no=0)
          l.b[i] = sum((delta.b-mu.b[i])*w.b)/sum(w.b)
          l.b[is.na(l.b)] = 0 # remove NAs
        }
        l = c(l.a,l.b)
        W.tilde = (X[,D.ind[d.ind]]-d0.c[j]) - l
        # uniform kernel
        Kh = (1/h[j])*ifelse(abs(X[,D.ind[d.ind]]-d0.c[j])/h[j]<=1, yes = 1/2, no = 0)
        W.hat = W.tilde - sum(W.tilde*Kh)/sum(Kh)
        Sn.hat = 1/n*sum(W.hat*(X[,D.ind[d.ind]]-d0.c[j])*Kh)
        est[j,d.ind] = 1/(n*Sn.hat)*sum(W.hat*R.hat*Kh)
        V.hat = sigma1.sq/(n^2*Sn.hat^2)*sum(W.hat^2*Kh^2)
        est.se[j,d.ind] = sqrt(V.hat)
      } # for each evaluation point

    } else {
      if (!treatment.SAM) {
        gamma = Initialization.step(X[,-D.ind[d.ind]], X[,D.ind[d.ind]], lambda = "CV.min", intercept = TRUE)$lasso.est
      } else {
        sam.D = cv.SAM(X[,-D.ind[d.ind]],X[,D.ind[d.ind]],lam.seq = lam.seq,quant.trans = T)
        D.hat = predict(sam.D$sam.final,newdata=X[,-D.ind[d.ind]])$values
      }
      for (j in 1:n.d0) {
        # calculate mu and delta
        if (!treatment.SAM) {
          mu = d0.c[j] - cbind(1,matrix(X[,-D.ind[d.ind]],n,p-1))%*%gamma
          delta= matrix(X[,D.ind[d.ind]],n,1) - cbind(1,matrix(X[,-D.ind[d.ind]],n,p-1))%*%gamma
        } else {
          mu = d0.c[j] - D.hat
          delta = matrix(X[,D.ind[d.ind]],n,1) - D.hat
        }
        # calculate l
        l = rep(0,n)
        for (i in 1:n) {
          # weight
          w = ifelse(abs(delta-mu[i])<=h[j],yes=1,no=0)
          l[i] = sum((delta-mu[i])*w)/sum(w)
          l[is.na(l)] = 0 # remove NAs
        }
        # center the decorrelation weights
        W.tilde = (X[,D.ind[d.ind]]-d0.c[j]) - l
        # uniform kernel
        Kh = (1/h[j])*ifelse(abs(X[,D.ind[d.ind]]-d0.c[j])/h[j]<=1, yes = 1/2, no = 0)
        W.hat = W.tilde - sum(W.tilde*Kh)/sum(Kh)
        Sn.hat = 1/n*sum(W.hat*(X[,D.ind[d.ind]]-d0.c[j])*Kh)
        est[j,d.ind] = 1/(n*Sn.hat)*sum(W.hat*R.hat*Kh)
        V.hat = sigma1.sq/(n^2*Sn.hat^2)*sum(W.hat^2*Kh^2)
        est.se[j,d.ind] = sqrt(V.hat)
      } # for each evaluation point
    }
  }

  # confidence interval
  # CI = rep(list(NA),2)
  # names(CI) = c("lower","upper")
  # CI[[1]] = est + qnorm(alpha/2)*est.se
  # CI[[2]] = est + qnorm(1-alpha/2)*est.se
  CI = rep(list(NA), length(D.ind))
  names(CI) = paste("f",as.character(D.ind),sep="")
  for(d.ind in 1:length(D.ind)){
    CI[[d.ind]] = matrix(NA, nrow=n.d0, ncol=2)
    colnames(CI[[d.ind]]) = c("lower", "upper")
    rownames(CI[[d.ind]]) = d0
    CI[[d.ind]][,1] = est[,d.ind] + qnorm(alpha/2)*est.se[,d.ind]
    CI[[d.ind]][,2] = est[,d.ind] + qnorm(1-alpha/2)*est.se[,d.ind]
  }
  returnList = list(est = est,
                    est.se = est.se,
                    CI = CI,
                    d0 = d0,
                    bw.save = bw.save,
                    sigma1.sq = sigma1.sq
  )
  returnList
}
