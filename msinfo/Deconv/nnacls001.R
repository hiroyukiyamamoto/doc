# non-negative alternating least squares
# Yamamoto et al.

library(osd)
library(nnls)

data(gcms1) # 71(クロマトの時間ポイント)×566(MSスペクトルのm/zポイント)
X <- gcms1

# number of component
com <- 3
maxiter <- 100

# initial value
C <- matrix(runif(nrow(X)*com),ncol=com)

E <- NULL
for (k in 1:maxiter){
  A <- NULL
  for(i in 1:ncol(X)){
    a <- nnls(C,X[,i])$x
    A <- cbind(A,a)
  }

  # normalized constraint
  for (i in 1:com){
    A[i,] <- A[i,]/as.numeric(sqrt(t(A[i,])%*%A[i,]))
  }
  
  C <- NULL
  for(i in 1:nrow(X)){
    c0 <- nnls(t(A),(X[i,]))$x
    C <- cbind(C,c0)
  }
  C <- t(C)
  
  # error
  E[k] <- norm(X-C%*%A);  
}

par(mfrow = c(3, 2)) 
plot(C[,1], type='l')
plot(A[1,], type='h')
plot(C[,2], type='l')
plot(A[2,], type='h')
plot(C[,3], type='l')
plot(A[3,], type='h')

#plot(E,type="l")
  