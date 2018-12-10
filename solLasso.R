#' Game: Calculate 4 number to 24.
#' @usage Calc24(c(a,b,c,d))
#' @param a A integer.
#' @param b A integer.
#' @param c A integer.
#' @param d A integer.
#' @return The final formula to calculate 24 using \code{a}, \code{b}, \code{c} and \code{d} by  +, -, *, and /.
#' @examples
#' Calc24(c(4,4,4,4))
#' [1] '((4*4)+4)+4=24'
#' #This is the formula to calculate 24.#
#'
#' Calc24(c(1,2,3,4))
#' [1] '((4*3)*2)*1=24'
#' #This is the formula to calculate 24.#
##'@export
l2n <- function(p){
  sqrt(t(p) %*% p)
}

##'@export
diffchange <- function(p, q){
  deltavec <- p - q
  sumvec <- p + q
  l2n(deltavec)/l2n(sumvec)
}

##'@export
softThreshold <- function(yy, lambda){
  sign(yy) * pmax(abs(yy) - lambda, 0)
}

##'@export
prox <- function(yy, lambda, t){
  softThreshold(yy=yy, lambda = lambda*t)
}

##'@export
admm <- function(b, rou, v, lambda){
  softThreshold(yy=b+1/rou*v, lambda = lambda/rou)
}
##'@export
scaled_admm <- function(b, rou, u, lambda){
  softThreshold(yy=b+u, lambda = lambda/rou)
}
##'@export
f <- function(b, y, x, lambda){1/2*l2n(y-x%*%b)^2 + lambda*sum(abs(b))}
##'@export
m <- function(b, y, x){1/2*l2n(y-x%*%b)^2}
##'@export
n <- function(b, y, x){-t(x)%*%(y-x%*%b)}
##'@export
g <- function(b, y, x, t, lambda){(b-prox(yy=b-t*n(b=b,y=y,x=x), lambda=lambda, t=t))/t}
##'@export
pgdLasso <- function(x,y,lambda,inib = NULL,t0 = 1,z = 1/2,error = 1e-16, silence = F){
  
  if(is.null(inib)){
    b <- replicate(ncol(x), 1)
  } else {
    b <- inib
  }
  b_old <- replicate(ncol(x),0)
  step <- 0

  while(diffchange(b,b_old)>error){
    t=t0
    b_old <- b
    
    gg <- g(b=b,y=y,x=x,t=t,lambda = lambda)
    l <- m(b=b-t*gg, x=x, y=y)
    r <- m(b=b,y=y,x=x)- t * t(n(b=b,y=y,x=x)) %*% gg + 1/2*t*l2n(gg)^2
    
    while(l > r) {
      t = t*z
      gg <- g(b=b,y=y,x=x,t=t,lambda = lambda)
      l <- m(b=b-t*gg, x=x, y=y)
      r <- m(b=b,y=y,x=x)- t * t(n(b=b,y=y,x=x)) %*% gg + 1/2*t*l2n(gg)^2
    }
    
    b <- softThreshold(yy=b+t*t(x)%*%(y-x%*%b),lambda = t*lambda)
    step <- step + 1
    if(silence==F){
      cat("Step:", step,"\n", "OFV:", f(b=b,x=x,y=y,lambda=lambda),"\n")
    }
  }
  list(NoofSteps = step, Beta = b)
}
##'@export
cdLasso <- function(x,y,lambda,inib = NULL,error = 1e-16,silence = F){
  
  if(is.null(inib)){
    b <- replicate(ncol(x), 1)
  } else {
    b <- inib
  }
  b_old <- replicate(ncol(x), 0)
  step <- 0
  
  while (diffchange(b,b_old)>error) {
    
    b_old <- b
    
    for (j in 1:ncol(x)) {
      
      xjty <- t(x[,j])%*%y
      xjtxipi <- t(x[,j])%*%rowSums(x[,-j]%*%b[-j])
      xjtxj <- t(x[,j])%*%x[,j]
      
      ifelse(xjty > xjtxipi + lambda, 
             b[j] <- (xjty-xjtxipi-lambda)/xjtxj, 
             ifelse(xjty < xjtxipi - lambda, 
                    b[j] <- (xjty - xjtxipi + lambda)/xjtxj, 
                    b[j] <- 0))
    }
    step <- step +1
    if(silence==F){
      cat("Step:", step,"\n", "OFV:", f(b=b,x=x,y=y,lambda=lambda),"\n")
    }
  }
  list(NoofSteps = step, Beta = b)
}
##'@export
admmLasso <- function(x,y,lambda,inib = NULL,error = 1e-16,rou = 1,silence = F,bthreshold = 1e-3){
  
  if(is.null(inib)){
    b <- replicate(ncol(x), 1)
  } else {
    b <- inib
  }
  b_old <- replicate(ncol(x),0)
  v <- replicate(ncol(x), 1)
  z <- replicate(ncol(x), 1)
  step <- 0
  xtx <- t(x)%*%x
  xtxr <- solve(t(x)%*%x + rou*diag(1, ncol(x), ncol(x)))
  xty <- t(x)%*%y
  
  while (diffchange(b,b_old)>error){
    
    b_old <- b
    b <- xtxr %*% (xty - v + rou*z)
    z <- admm(b=b, rou=rou, v=v, lambda=lambda)
    v <- v + rou*(b - z)
    step = step + 1
    if(silence==F){
      cat("Step:", step,"\n", "OFV:", f(b=b,x=x,y=y,lambda=lambda),"\n")
    }
  } 
  b[b < bthreshold] <- 0
  list(NoofSteps = step, Beta = b)
}
##'@export
scadmmLasso <- function(x,y,lambda,inib = NULL,error = 1e-16,rou = 1,silence = F,bthreshold = 1e-3){
  
  if(is.null(inib)){
    b <- replicate(ncol(x), 1)
  } else {
    b <- inib
  }
  b_old <- replicate(ncol(x),0)
  u <- replicate(ncol(x), 1)
  z <- replicate(ncol(x), 1)
  step <- 0
  xtx <- t(x)%*%x
  xtxr <- solve(t(x)%*%x + rou*diag(1, ncol(x), ncol(x)))
  xty <- t(x)%*%y
  
  while (diffchange(b,b_old)>error){
    
    b_old <- b
    b <- xtxr %*% (xty - rou*u + rou*z)
    z <- scaled_admm(b=b, rou=rou, u=u, lambda=lambda)
    u <- u + b - z
    step = step + 1
    if(silence==F){
      cat("Step:", step,"\n", "OFV:", f(b=b,x=x,y=y,lambda=lambda),"\n")
    }
  } 
  b[b < bthreshold] <- 0
  list(NoofSteps = step, Beta = b)
}




library(ElemStatLearn)
data("prostate", package = "ElemStatLearn")
summary(prostate)

x <- replicate(ncol(prostate)-1, c(1:nrow(prostate)))
x[,2:9] <- as.matrix(prostate[,1:8])
x[,2:9] <- scale(as.matrix(prostate[,1:8]))
x[,1] <- 1
y <- replicate(1, c(1:nrow(prostate)))
y <- as.matrix(prostate[,9])

pgdLasso(x=x,y=y,lambda = 10)
cdLasso(x=x,y=y,lambda = 10)
admmLasso(x=x,y=y,lambda = 10, silence = T)
scadmmLasso(x=x,y=y,lambda = 10, silence = T)
