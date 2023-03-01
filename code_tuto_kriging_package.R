# Code for the tutorial seminar on Kriging and Bayesian Optimization
# Package DiceKriging
#
# Author : Tanguy Appriou
# Date : 2 March 2023

library(DiceDesign)
library(DiceKriging)
library(DiceEval)
library(DiceOptim)


### Branin test function ###


branin <- function(X) {
  if (!is.matrix(X)) X <- matrix(X,ncol=2) #Transform X into a matrix if it is a vector
  if (ncol(X) != 2) {
    if (nrow(X) != 2) stop(paste("The dimension must be ", 2))
    else X <- t(X)
  }
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

# Construct and plot the design plan (LHS)
d <- 2 #Dimension
n_train <- 10 #Number of training points
set.seed(123)
X_train_LHS <- lhsDesign(n=n_train, dimension=d, randomized=FALSE, seed=123)$design
X_train <- maximinESE_LHS(X_train_LHS)$design #Optimize the LHS
x11()
plot(X_train, pch=4, lwd=3, cex=1.5, col="black", 
     xlim=0:1, ylim=0:1, xlab=expression(x[1]), ylab=expression(x[2]), main="Initial design plan (Optimized LHS)")

# Function values at the training points
Y_train <- branin(X_train)



### Optimization of hyperparameters and building the Kriging model ###

model <- km(design=X_train, response=Y_train, covtype="gauss", estim.method="MLE", optim.method="BFGS",
            lower=rep(10^-1,d), upper=rep(10,d), nugget=NULL)



### Predictions ###

pred <- predict(object=model, newdata=X_grid, type="UK", se.compute=TRUE, checkNames=FALSE)
Y_pred <- pred$mean
var_pred <- (pred$sd)^2

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



### EGO ###

n_steps <- 15 #Number of EGO steps
EGO <- EGO.nsteps(model=model, fun=branin, nsteps=n_steps, lower=rep(0,d), upper=rep(1,d), control=list(pop.size=20, max.generations=50, trace=0))

# Predictions after the EGO iterations
pred_afterEGO <- predict(object=EGO$lastmodel, newdata=X_grid, type="UK", se.compute=TRUE, checkNames=FALSE)
Y_pred_afterEGO <- pred_afterEGO$mean

# Plot
x11()
filled.contour(matrix(Y_pred_afterEGO, neval, neval), xlim=0:1, ylim=0:1, nlevels=30, color.palette = palette, 
               xlab=expression(x[1]), ylab=expression(x[2]), main="Mean prediction",
               plot.axes={axis(1); axis(2); points(X_train, pch=4, cex=1.5, lwd=3, col="black");
                 points(EGO$par, pch=4, cex=1.5, lwd=3, col="red");
                 legend("topright", legend=c("Initial Design points","Additional points"), pch=rep(4,2), lty=rep(NA,2), pt.cex=rep(1.5,2), lwd=rep(3,2), col=c("black","red"))})
