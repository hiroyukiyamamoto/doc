# Regularized alternating least squares
# Yamamoto et al.

library(osd)

data(gcms1) # 71(クロマトの時間ポイント)×566(MSスペクトルのm/zポイント)
X <- gcms1

# number of component
com <- 3
lambda <- 0
maxiter <- 1000

# initial value
C <- matrix(runif(nrow(X)*com),ncol=com)

E <- NULL
for (k in 1:maxiter){
  A <- solve(t(C)%*%C+lambda*diag(1,com))%*%t(C)%*%X
  A[A<0] <- 0
                                      
  # normalized constraint
  for (i in 1:com){
    A[i,] <- A[i,]/as.numeric(sqrt(t(A[i,])%*%A[i,]))
  }
      
  C <- X%*%t(A)%*%solve(A%*%t(A)+lambda*diag(1,com))
  C[C<0] <- 0
  
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

  