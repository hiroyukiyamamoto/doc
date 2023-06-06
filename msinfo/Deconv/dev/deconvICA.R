### ICA
library(ica)


deconvICA <- function(Y0){
  
  com <- 5
  icacom <- ica(t(Y0),com)
  
  C <- matrix(NA,com,nrow(Y))
  for(k in 1:com){
    S <- icacom$M[,k]
    R <- cor(Y0,S)
    
    if (R[which.max(abs(R))]<0){
      S <- -S
    }
    S[S<0] <- 0
    C[k,] <- S
  }

  # number of component
  lambda <- 0.01
  maxiter <- 5
  
  # initial value
  C <- t(C)
  
  X <- Y
  
  E <- NULL
  for (k in 1:maxiter){
    A <- solve(t(C)%*%C+lambda*diag(1,com))%*%t(C)%*%X
    A[A<0] <- 0
    
    # normalized constraint
    #for (i in 1:com){
    #  A[i,] <- A[i,]/as.numeric(sqrt(t(A[i,])%*%A[i,]))
    #}
    
    C <- X%*%t(A)%*%solve(A%*%t(A)+lambda*diag(1,com))
    C[C<0] <- 0
    
    # error
    E[k] <- norm(X-C%*%A);
  }
  
  A_fin <- NULL
  for(i in 1:com){
    if(sum(A[i,]!=0)>=1){
      A_fin <- rbind(A_fin, A[i,])  
    }
  }
  
  return(A_fin)
}
