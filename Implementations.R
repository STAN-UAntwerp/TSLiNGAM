library(pracma); library(dccpp); library(glmnet)
# DirectLiNGAM code based on Matlab code accompanying the original paper of Shimizu et al. (2011)

DirectLingam <- function(X, M = NULL, indep_meas = 1, scale = 0) {
  # X : n x p input data matrix
  # M : p x p matrix holding prior knowledge on structure
  
  # indep_meas = 1 : kernel based independence measure
  # indep_meas = 2 : distance correlation
  
  Xorig <- X
  n <- nrow(X)
  p <- ncol(X)
  
  # centering (and scaling) by subtracting the column means
  if (scale == 1){
    X <- scale(X) 
  }else{
    X <- scale(X, scale = FALSE) 
  }
  
  
  if (is.null(M)) {
    M <- matrix(-1, p, p) # no prior knowledge
  }
  
  # erase diagonal elements of M
  diag(M) <- 0
  
  ##########
  # STEP 1 #
  ##########
  
  K <- rep(0, p) # ordered list of variables
  U_K <- 1:p # variables left (U\K)
  
  ##########
  # STEP 2 #
  ##########
  
  for (m in 1:(p-1)) {
    # find exogenous by using M
    exogenous = which(rowSums(M == 0, na.rm = TRUE) == p - m +1)  # length(U_K) == p-m+1
    if (length(exogenous) == 0) {
      endogenous = which(rowSums(M == 1) > 0)
      candidates = setdiff(U_K, endogenous)
    } else {
      candidates = exogenous
    }
    
    # calculate residuals R
    R = computeR(X, candidates, U_K, M)
    
    # skip exogenous finding if it is found
    if (length(exogenous) > 0) { 
      index = exogenous[1]
    } else if (length(candidates) == 1){
      index = candidates
    }
    else {
      # find exogenous
      index = findindex(X, R, candidates, U_K, indep_meas)
    }
    
    K[m] = index
    U_K = U_K[-which(U_K == index)]
    M[,index] = NA
    M[index,] = NA;
    
    X = t(R[, , index]);
  }
  
  ##########
  # STEP 3 #
  ##########
  
  K[p] = U_K;
  
  ##########
  # STEP 4 #
  ##########
  
  ols.out <- ols(Xorig, K) # perform final regression to estimate strengths and constants
  
  return(list(B = ols.out$Bols,
              stde = ols.out$olsdisturbancestd,
              ci = ols.out$cols,
              K = K))
}


TSLiNGAM <- function(X, M = NULL, indep_meas = 1, slope = 1, scale = 0) {
  # X : n x p input data matrix
  # M : p x p matrix holding prior knowledge on structure
  
  # indep_meas = 1 : kernel based independence measure
  # indep_meas = 2 : distance correlation
  
  # slope = 1 : Theil-Sen
  # slope = 2 : Repeated Median
  
  # for the TSLiNGAM method: indep_meas = 1 and slope = 1
  
  Xorig <- X
  n <- nrow(X)
  p <- ncol(X)
  
  # centering: use robust location estimate instead of mean
  if (scale == 1){
    s <- apply(X, 2, Qn)
    s[s==0]=1 # if Qn is 0, make it 1 so no scaling is done on that variable
    X <- scale(X, center = robustbase::colMedians(X), scale = s) 
  }else{
    X <- scale(X, center = robustbase::colMedians(X), scale = FALSE)
  }
  
  if (is.null(M)) {
    M <- matrix(-1, p, p)
  }
  
  # erase diagonal elements of M
  diag(M) <- 0
  
  ##########
  # STEP 1 #
  ##########
  
  K <- rep(0, p)
  U_K <- 1:p
  
  ##########
  # STEP 2 #
  ##########
  
  for (m in 1:(p-1)) {
    # find exogenous by using M
    exogenous = which(rowSums(M == 0) == p - m + 1) # length(U_K) == p-m+1
    if (length(exogenous) == 0) {
      endogenous = which(rowSums(M == 1) > 0)
      candidates = setdiff(U_K, endogenous)
    } else {
      candidates = exogenous
    }
    
    # calculate residuals R
    R = computeR_rob(X, slope, candidates, U_K, M)
    
    if (length(candidates) == 1) {
      index = candidates
    } else {
      if (indep_meas == 1){
        index = findindex( X, R, candidates, U_K, indep_meas)
      } else if (indep_meas == 2){ 
        index = findindex( X, R, candidates, U_K, indep_meas)
      } else if (indep_meas == 3){
        index = findindex( X, R, candidates, U_K, indep_meas)
      } else if(indep_meas == 4){
        index = findindex( X, R, candidates, U_K, indep_meas) 
      }
    }
    K[m] = index
    U_K = U_K[-which(U_K == index)]
    M[,index] = NA
    M[index,] = NA;
    
    X = t(R[, , index]);
  }
  
  ##########
  # STEP 3 #
  ##########
  
  K[p] = U_K;
  
  ##########
  # STEP 4 #
  ##########
  
  ols.out <- ols(Xorig, K) # perform final regression to estimate strengths and constants
  
  return(list(B = ols.out$Bols,
              stde = ols.out$olsdisturbancestd,
              ci = ols.out$cols,
              K = K))
}


ols <- function(X, K) {
  # K the estimated causal order of the variables
  # X the observed n x p data matrix
  
  p <- ncol(X)
  
  # permute the variables to the causal order
  X = X[, K];
  
  # subtract out the mean
  Xm <- colMeans(X)
  X <- scale(X, center = Xm, scale = FALSE)
  
  # calculate covariance matrix
  C = cov(X)
  
  qlDecomp = function(A) {
    # only for square matrices
    # Finds Q and L s.t. QL = A with 
    # Q orthonormal columns
    # and L a lower triangular matrix
    p = nrow(A)
    iI = apply(diag(p), 2, rev)
    qr.out = qr(iI %*% A %*% iI)
    Q = iI %*%qr.Q(qr.out)%*% iI
    L = iI %*%qr.R(qr.out)%*% iI
    return(list(Q = Q, L = L))
  }
  
  # do QL decomposition on the inverse square root of C
  ei <- eigen(C)
  V <- ei$vectors
  Csqrtinv <- V %*% diag(1 / sqrt(ei$values)) %*% t(V)
  QL.out   <- qlDecomp(Csqrtinv)
  L        <- QL.out$L
  
  # the estimated disturbance-stds are one over the abs of the diag of L
  olsdisturbancestd = 1 / diag(abs(L))
  
  # normalize rows of L to unit diagonal
  L = t(scale(t(L), FALSE, diag(L)))
  
  # calculate corresponding B
  Bols = diag(p) - L
  
  # also calculate constants: corresponds to the e-vector
  cols = L %*% Xm
  
  # permute back to original variable order
  Bols[K, K] = Bols
  olsdisturbancestd[K] = olsdisturbancestd
  cols[K] = cols;
  
  return(list(Bols = Bols, olsdisturbancestd = olsdisturbancestd, cols = cols))
}


findindex = function(X, R, candidates, U_K, indep_meas = 1) {
  
  p <- ncol(X)
  T <- matrix(NA, 1, p)
  minT <- -1
  
  for (j in candidates) {
    if (minT == -1){
      T[j] = 0
      for (i in setdiff(U_K, j)){
        if (indep_meas == 1){
          J = kernel_based_independence_measure(cbind(R[i,,j],X[,j]))
        }else if(indep_meas == 2){
          J = dccpp::dcor(R[i,,j],X[,j])
        }
        T[j] = T[j] + J
      }
      minT = T[j]
    } else {
      T[j] = 0
      for (i in setdiff(U_K, j)){
        if (indep_meas == 1){
          J = kernel_based_independence_measure(cbind(R[i,,j],X[,j]))
        }else if (indep_meas == 2){
          J = dccpp::dcor(R[i,,j],X[,j])
        }
        T[j] = T[j] + J
        if (T[j] > minT){
          T[j] = Inf
          break
        }
      }
      minT <- min(c(T[j],minT))
    }
  }
  dummy <- which(T == (min(T, na.rm = TRUE)))
  return(dummy[1])
}


kernel_based_independence_measure = function(X){
  # Code based on Matlab code by Francis R. Bach, 2002.
  
  n <- nrow(X)
  p <- ncol(X)
  
  X <- t(X) # Matlab convention
  
  if (n <= 1000){
    sigma = 1
    kappa = 2e-2
  }
  else {
    sigma = 1/2
    kappa = 2e-3
  }
  eta = kappa*(1e-2)
  
  Us <- matrix(list(), 1, p)
  Lambdas <- matrix(list(), 1, p)
  Drs <- matrix(list(), 1, p)
  sizes <- vector()
  
  for(i in 1:p){
    chol <- chol_Gauss(X[i,]/sigma,1,n*eta)
    G <- chol$G
    pvec <- chol$pvec
    s <- sort(pvec, decreasing = FALSE, index.return = TRUE)
    pvec <- s$ix
    G <- scale(G[pvec,], scale = FALSE)
    eigens <- eigen(t(G) %*% G)
    A <- eigens$vectors
    D <- eigens$values
    indexes = which(D >= n*eta & Im(D)==0)
    s2 <- sort(D[indexes],decreasing = FALSE, index.return = TRUE)
    if(length(D[indexes]) != 1){
      order <- s2$ix
      order <- rev(order)
    }else{
      order <- 1
    }
    neig <- length(indexes)
    if (isempty(indexes)){
      indexes = c(1)
    } else {
      indexes <- indexes[order[1:neig]]
    }
    D <- D[indexes]
    if(length(D[indexes]) != 1){
      V <- G %*% (A[,indexes] %*% diag(sqrt(1/D)))
    }else{
      V <- G %*% (A[,indexes] * sqrt(1/D))
    }
    
    Us[[1,i]] <- V
    Lambdas[[1,i]] <- D
    Dr <- D
    for (j in 1:length(D)){
      Dr[j] <- D[j]/(n*kappa + D[j])
    }
    Drs[[1,i]] <- Dr
    if(length(D[indexes]) != 1){
      sizes[i] <- length(Drs[[1,i]])
    }else{
      sizes[i] <- 1
    }
  }
  Rkappa <- diag(sum(sizes))
  starts <- cumsum(c(1,sizes))
  starts <- starts[1:p]
  for ( i in 2:p){
    for ( j in 1:(i-1) ){
      if(length(Drs[[1,i]])== 1 & length(Drs[[1,j]]) ==1){
        newbottom <- Drs[[1,i]] %*% (t(Us[[1,i]]) %*% Us[[1,j]]) %*% Drs[[1,j]]
      }else if (length(Drs[[1,i]])== 1){
        newbottom <- Drs[[1,i]] %*% (t(Us[[1,i]]) %*% Us[[1,j]]) %*% diag(Drs[[1,j]])
      }else if (length(Drs[[1,j]])== 1){
        newbottom <- diag(Drs[[1,i]]) %*% (t(Us[[1,i]]) %*% Us[[1,j]]) %*% Drs[[1,j]]
      }else{
        newbottom <- diag(Drs[[1,i]]) %*% (t(Us[[1,i]]) %*% Us[[1,j]]) %*% diag(Drs[[1,j]])
      }
      Rkappa[starts[i]:(starts[i]+sizes[i]-1),starts[j]:(starts[j]+sizes[j]-1)] = newbottom
      Rkappa[starts[j]:(starts[j]+sizes[j]-1),starts[i]:(starts[i]+sizes[i]-1)] = t(newbottom)
    }
  }
  D = det(Rkappa)
  J = -0.5*log(D)
  
  return(J)
}


chol_Gauss <- function(X,sigma,tol){
  # Code based on Matlab code by Francis R. Bach, 2002.
  
  X = t(X)
  n <- ncol(X) # Matlab convention
  pvec <- 1:n
  diagG <- rep(1,n)
  i <- 1
  G <- matrix(, nrow = n, ncol = 0)
  
  while( (sum(diagG[i:n]) > tol) & (i <= n) ){
    
    G <- cbind(G,rep(0,n))
    
    if (i > 1){
      jast <- which(diagG[i:n] == (max(diagG[i:n], na.rm = TRUE)))
      jast <- jast[1]
      jast = jast + i - 1
      pvec[c(i,jast)] = pvec[c(jast,i)]
      G[c(i,jast),1:i] = G[c(jast,i),1:i]
    }
    else{
      jast = 1
    }
    
    G[i,i] = sqrt(diagG[jast])
    
    sqdist = function(a,b){
      if (pracma::size(a)[1] == 1){ # for a vector
        aa = a*a
        a = matrix(a,byrow=T)
      }else{ # for a matrix
        aa = colSums(a*a)
      }
      if (pracma::size(b)[1] == 1){
        bb = b*b
        b = matrix(b,nrow=1)
      }else{
        bb = colSums(b*b)
      }
      ab = a %*% b
      d = abs(repmat(matrix(aa,byrow=T),1,length(bb)) + repmat(bb,length(aa),1) - 2*ab)
      return(d)
    }
    
    if (i < n){
      newAcol = exp( -(0.5/sigma^2)*sqdist(X[,pvec[(i+1):n]],X[,pvec[i]]) )
      if (i > 2){
        G[(i+1):n,i] = 1/G[i,i] * (newAcol - G[(i+1):n,1:(i-1),drop=FALSE] %*% t(G[i,1:(i-1),drop=FALSE]))
      } else if (i == 2){
        G[3:n,2] = 1/G[2,2] * (newAcol - G[3:n,1,drop=FALSE]*G[2,1])
      }else{ # i=1
        G[2:n,1] = 1/G[1,1] * newAcol
      }
    }
    if (i < n){
      if (i == 1){
        diagG[(i+1):n] = matrix(rep(1,n-i),nrow=n-i) - matrix(G[2:n,1]^2,nrow=n-1)
      }else if (i == n-1) {
        diagG[n] = 1 - sum(G[n,1:i]^2)
      }
      else{
        diagG[(i+1):n] = matrix(rep(1,n-i),nrow=n-i) - matrix(rowSums(G[(i+1):n,1:i]^2),nrow=n-i)
      }
    }
    i = i + 1
  }
  return(list(G = G, pvec = pvec))
}


computeR <- function(X, candidates, U_K, M) {
  # computes the residuals
  
  n <- nrow(X)
  p <- ncol(X)
  R <- array(0, dim = c(p, n, p))
  
  Cov <- cov(X)
  
  for (j in candidates) {
    for (i in setdiff(U_K, j)) {
      # skip residue calculation by using M
      if (M[i,j] == 0) {
        R[i,,j] = X[ ,i]
      } else {
        R[i,,j] = X[, i] - (Cov[i, j] / Cov[j,j]) * X[,j]
      }
    }
  }
  return(R)
}


computeR_rob <- function(X, slope, candidates, U_K, M) {
  # computes the residuals using the Theil-Sen slope or the Repeated Median
  
  n <- nrow(X)
  p <- ncol(X)
  R <- array(0, dim = c(p, n, p))
  
  for (j in candidates) {
    for (i in setdiff(U_K, j)) {
      # skip residue calculation by using M
      if (M[i,j] == 0) {
        R[i, ,j] = X[ ,i]
      } else {
        if (slope == 1){
          robslope = robslopes::TheilSen(X[, j], X[, i], verbose = FALSE)
        }
        else if (slope == 2){
          robslope = robslopes::RepeatedMedian(X[, j], X[, i], verbose = FALSE)
        }
        R[i,,j] = X[, i] - robslope$slope * X[,j]
      }
    }
  }
  return(R)
}


Bprune_adlasso <- function(X, k, nCV = 5, B){
  # Prunes the adjacency matrix B to remove redundant edges
  
  # X : data matrix
  # k : estimated causal order
  # nCV: nCV-fold cross-validation, default = 5
  # B : non pruned B matrix found with OLS
  
  p <- ncol(X)
  
  Xp <- X[,k] # permute to causal order
  Bp <- B[k,k] # permute to causal order
  
  Badap <- zeros(dim(B)[1])
  
  # for first two variables we can't do adaptive lasso with glmnet (multiple predictors needed)
  # hence we just look whether the regression with 1 predictor achieves a minimal R^2
  data12 <- data.frame(Xp[,1:2])
  colnames(data12) <- c("var1","var2")
  model12 <- lm(var2~var1 - 1,data = data12)
  rsquared <- summary(model12)$r.squared
  if (rsquared > 0.05){ # minimal requirement for model
    Badap[2,1] = Bp[2,1]
  }
  
  # for other variables: adaptive lasso
  if (p > 2){
    for (i in 3:p){
      resp = Xp[,i]
      pred = Xp[,1:(i-1)]
      OLS_zelf = Bp[i,1:(i-1)]
      
      # now adaptive lasso with 3 gamma's: 0.5, 1, 2
      alasso1_cv_05 <- cv.glmnet(x = pred, y = resp, alpha = 1, intercept=F,
                                 penalty.factor = 1 / (abs(OLS_zelf))^0.5)
      alasso1_cv_1  <- cv.glmnet(x = pred, y = resp, alpha = 1, intercept=F,
                                 penalty.factor = 1 / (abs(OLS_zelf))^1)
      alasso1_cv_2  <- cv.glmnet(x = pred, y = resp, alpha = 1, intercept=F,
                                 penalty.factor = 1 / (abs(OLS_zelf))^2)
      s_05 <- alasso1_cv_05$lambda.min
      s_1  <- alasso1_cv_1$lambda.min
      s_2  <- alasso1_cv_2$lambda.min
      best_alasso_coef_05 <- coef.glmnet(alasso1_cv_05, s = s_05)
      best_alasso_coef_1  <- coef.glmnet(alasso1_cv_1,  s = s_1)
      best_alasso_coef_2  <- coef.glmnet(alasso1_cv_2,  s = s_2)
      res_05 <- sum((resp - matrix(best_alasso_coef_05@x,nrow=1)%*%t(pred[,best_alasso_coef_05@i]))^2)
      res_1  <- sum((resp - matrix(best_alasso_coef_1@x, nrow=1)%*%t(pred[,best_alasso_coef_1@i ]))^2)
      res_2  <- sum((resp - matrix(best_alasso_coef_2@x, nrow=1)%*%t(pred[,best_alasso_coef_2@i ]))^2)
      min <- which.min(c(res_05, res_1, res_2)) # choose best gamma
      if(min == 1){
        Badap[i,best_alasso_coef_05@i] <- Bp[i,best_alasso_coef_05@i]
      }else if (min == 2){
        Badap[i,best_alasso_coef_1@i]  <- Bp[i,best_alasso_coef_1@i]
      }else if (min == 3){
        Badap[i,best_alasso_coef_2@i]  <- Bp[i,best_alasso_coef_2@i]
      }
    }
  }
  
  # Badap inverse permutation
  inverse_perm = sort(k,index.return = T)$ix
  Badap = Badap[inverse_perm,inverse_perm]
  
  return(Badap)
}

