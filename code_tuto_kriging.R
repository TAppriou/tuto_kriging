# Code for the tutorial seminar on Kriging and Bayesian Optimization
#
# Author : Tanguy Appriou
# Date : 2 March 2023


### Branin test function ###


branin <- function(X) {
  if (!is.matrix(X)) X <- matrix(X,ncol=1) #Transform X into a matrix if it is a vector
  if (ncol(X) != 2) stop(paste("X must have", 2, " columns"))
  x1 <- X[,1]*15 - 5
  x2 <- X[,2]*15
  res <- (x2 - 5.1/(4*pi^2) * (x1^2) + 5/pi*x1 - 6)^2 + 10*(1 - 1/(8*pi))*cos(x1) + 10 + 5*(x1+5)/15
  return(res)
}


# Plot the true Branin function on a grid
neval <- 201 #Number of points for plot is neval^d
X_grid <- as.matrix(expand.grid(seq(0,1,length.out=neval), seq(0,1,length.out=neval)))
Y_true <- branin(X_grid)
min_glob_branin <- c(-3.1762, 12.3582) #Global minimum of the Branin function
min_glob_branin <- (min_glob_branin + c(5,0))/15 #Normalization of the minima
min_loc_branin <- matrix(c(3.1069, 2.3023, 9.3901, 2.4461), nrow=2, ncol=2, byrow=TRUE) #Local minima of the Branin function
min_loc_branin <- (min_loc_branin + matrix(c(5,0),2,2,byrow=TRUE))/15 #Normalization of the minima

palette <- colorRampPalette(c("#0000A2", "#0017B9", "#002ED1", "#0046E8", "#005DFF", "#0074FF", "#008BFF", "#00A2FF", "#17B9FF", "#2ED1FF", "#46E8FF", "#5DFFFF", "#74FFE8", "#8BFFD1", "#A2FFB9", "#B9FFA2", "#D1FF8B", "#E8FF74", "#FFFF5D", "#FFE846", "#FFD12E", "#FFB917", "#FFA200", "#FF8B00", "#FF7400", "#FF5D00", "#E84600", "#D12E00", "#B91700", "#A20000"))
x11()
filled.contour(matrix(Y_true, neval, neval), xlim=0:1, ylim=0:1, nlevels=30, color.palette = palette, 
               xlab=expression(x[1]), ylab=expression(x[2]), main="True Branin function",
               plot.axes={axis(1); axis(2); points(t(min_glob_branin), pch=4, cex=1.5, lwd=3, col="red");
                 points(min_loc_branin, pch=4, cex=1.5, lwd=3, col="black");
                 legend("topright", legend=c("Global minimum","Local minima"), pch=rep(4,2), lty=rep(NA,2), pt.cex=rep(1.5,2), lwd=rep(3,2), col=c("red","black"))})



### Initial Design plan ###

# Function for a LHS design plan
lhs_design <- function(d,n_train) {
  design <- matrix(0, nrow=n_train, ncol=d)
  for (i in 1:d) {
    permut <- sample(1:n_train)
    design[,i] <- (permut - 0.5)/n_train
  }
  return(design)
}

# Construct and plot the design plan (LHS)
d <- 2 #Dimension
n_train <- 10 #Number of training points
set.seed(123) #Random seed
X_train_LHS <- lhs_design(d,n_train)
x11()
plot(X_train_LHS, pch=4, lwd=3, cex=1.5, col="black", 
     xlim=0:1, ylim=0:1, xlab=expression(x[1]), ylab=expression(x[2]), main="Initial design plan (LHS)")

# Optimize the LHS
library(DiceDesign)
X_train <- maximinESE_LHS(X_train_LHS)$design
x11()
plot(X_train, pch=4, lwd=3, cex=1.5, col="black", 
     xlim=0:1, ylim=0:1, xlab=expression(x[1]), ylab=expression(x[2]), main="Initial design plan (Optimized LHS)")

# Function values at the training points
Y_train <- branin(X_train)



### Optimization of hyperparameters ###

# Gaussian covariance function
cov_fun <- function(x,y,theta) {
  return( exp(-0.5 * sum( (x-y)^2/theta^2 )) )
}

# Create the covariance matrix between inputs
covmatrix_fun <- function(X,theta) {
  n <- dim(X)[1]
  K <- matrix(0,nrow=n,ncol=n)
  for (i in 1:n) {
    for (j in 1:i) {
      K[i,j] <- cov_fun(X[i,],X[j,],theta)
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}

# Covariance between new points and inputs
crosscovmatrix_fun <- function(X_train,theta,X_test) {
  n_train <- dim(X_train)[1]
  n_test <- dim(X_test)[1]
  K <- matrix(0,nrow=n_test,ncol=n_train)
  for (i in 1:n_test) {
    for (j in 1:n_train) {
      K[i,j] <- cov_fun(X_test[i,],X_train[j,],theta)
    }
  }
  return(K)
}

# pth partial derivative of the covariance matrix
dcovmatrix_fun <- function(X,theta,p) {
  n <- dim(X)[1]
  dK <- matrix(0,nrow=n,ncol=n)
  for (i in 1:n) {
    for (j in 1:i) {
      dK[i,j] <- (X[i,p]-X[j,p])^2 / theta[p]^3 * cov_fun(X[i,],X[j,],theta)
      dK[j,i] <- dK[i,j]
    }
  }
  return(dK)
}


# Function to optimize the Kriging length-scale hyperparameters
theta_krg_MLE <- function(X_train,Y_train) {
  
  n_train <- dim(X_train)[1]
  d <- dim(X_train)[2]
  
  ConcLogLik <- function(theta) {
    
    CholK <- chol.default(covmatrix_fun(X_train,theta))
    CholK_1 <- forwardsolve(t(CholK),matrix(1,nrow=n_train,ncol=1))
    mu <- crossprod(CholK_1, forwardsolve(t(CholK),matrix(Y_train,nrow=n_train,ncol=1))) / crossprod(CholK_1)
    CholK_Y <- forwardsolve(t(CholK),matrix(Y_train-c(mu),nrow=n_train,ncol=1))
    sigma2 <- 1/n_train * crossprod(CholK_Y)
    
    Loglik <- n_train * log(2*pi) + 2*sum(log(diag(CholK))) + n_train*log(sigma2) + n_train
    return(Loglik)
  }
  
  gradConcLogLik <- function(theta) {
    
    CholK <- chol.default(covmatrix_fun(X_train,theta))
    CholK_1 <- forwardsolve(t(CholK),matrix(1,nrow=n_train,ncol=1))
    mu <- crossprod(CholK_1, forwardsolve(t(CholK),matrix(Y_train,nrow=n_train,ncol=1))) / crossprod(CholK_1)
    CholK_Y <- forwardsolve(t(CholK),matrix(Y_train-c(mu),nrow=n_train,ncol=1))
    Kinv_Y <- backsolve(CholK,CholK_Y)
    sigma2 <- 1/n_train * crossprod(CholK_Y)
    
    grad <- rep(0,length=d)
    for (p in 1:d) {
      dK <- dcovmatrix_fun(X_train,theta,p)
      t1 <- - t(Kinv_Y) %*% dK %*% Kinv_Y / sigma2
      t2 <- sum(diag(chol2inv(CholK) %*% dK))
      grad[p] <- t1 + t2
    }
    return(grad)
  }
  
  lower_bound <- 10^-1
  upper_bound <- 10
  theta_init <- runif(2, min=lower_bound, max=upper_bound)
  opt <- optim(par=theta_init, ConcLogLik, gradConcLogLik, method="L-BFGS-B", lower=lower_bound, upper=upper_bound)
  
  return(opt$par)
}


# Optimize the hyperparameters
theta_krg <- theta_krg_MLE(X_train,Y_train)



### Kriging predictions ###


# Function for the Kriging predictions
pred_krg <- function(X_train,Y_train,theta,X_test) {
  
  n_train <- dim(X_train)[1]
  n_test <- dim(X_test)[1]
  Y_pred <- rep(0,length=n_test) #Mean
  var_pred <- rep(0,length=n_test) #Variance
  
  CholK <- chol.default(covmatrix_fun(X_train,theta))
  CholK_1 <- forwardsolve(t(CholK),matrix(1,nrow=n_train,ncol=1))
  mu <- crossprod(CholK_1, forwardsolve(t(CholK),matrix(Y_train,nrow=n_train,ncol=1))) / crossprod(CholK_1)
  CholK_Y <- forwardsolve(t(CholK),matrix(Y_train-c(mu),nrow=n_train,ncol=1))
  sigma2 <- 1/n_train * crossprod(CholK_Y)
  K_Xtest_X <- crosscovmatrix_fun(X_train,theta,X_test)
  CholK_K <- forwardsolve(t(CholK), t(K_Xtest_X))
  
  # Mean
  Y_pred <- matrix(mu, nrow=n_test, ncol=1) + crossprod(CholK_K, CholK_Y)
  # Variance
  var_pred <- c(sigma2) * (rep(1,length=n_test) - colSums(CholK_K^2))
  
  res <- data.frame(mean=Y_pred, var=var_pred)
  return(res)
}

# Predictions on the grid
pred <- pred_krg(X_train,Y_train,theta_krg,X_grid)
Y_pred <- pred$mean
var_pred <- pred$var

# Plot
x11()
filled.contour(matrix(Y_pred, neval, neval), xlim=0:1, ylim=0:1, nlevels=30, color.palette = palette, 
               xlab=expression(x[1]), ylab=expression(x[2]), main="Mean prediction",
               plot.axes={axis(1); axis(2); points(X_train, pch=4, cex=1.5, lwd=3, col="black");
                 legend("topright", legend=c("Design points"), pch=rep(4,1), lty=rep(NA,1), pt.cex=rep(1.5,1), lwd=rep(3,1), col=c("black"))})

x11()
filled.contour(matrix(var_pred, neval, neval), xlim=0:1, ylim=0:1, nlevels=30, color.palette = palette, 
               xlab=expression(x[1]), ylab=expression(x[2]), main="Prediction variance",
               plot.axes={axis(1); axis(2); points(X_train, pch=4, cex=1.5, lwd=3, col="black");
                 legend("topright", legend=c("Design points"), pch=rep(4,1), lty=rep(NA,1), pt.cex=rep(1.5,1), lwd=rep(3,1), col=c("black"))})

