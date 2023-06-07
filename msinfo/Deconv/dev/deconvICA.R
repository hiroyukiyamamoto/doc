### ICA
library(ica)

source("C:/R/function/matchedfilter.R")


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
  
  E <- NULL;A_ALL <- NULL;C_ALL <- NULL
  for (k in 1:maxiter){
    A <- solve(t(C)%*%C+lambda*diag(1,com))%*%t(C)%*%X
    A[A<0] <- 0
    A_ALL[[k]] <- A
    
    # normalized constraint
    #for (i in 1:com){
    #  A[i,] <- A[i,]/as.numeric(sqrt(t(A[i,])%*%A[i,]))
    #}
    
    C <- X%*%t(A)%*%solve(A%*%t(A)+lambda*diag(1,com))
    C[C<0] <- 0
    C_ALL[[k]] <- C
    
    # error
    E[k] <- norm(X-C%*%A);
  }
  
  # errorが最小のmaxiterを取得して、Aとする
  A <- A_ALL[[which.min(E)]]
  C <- C_ALL[[which.min(E)]]
  
  ## (単純にALSを使う場合は、ICAの段階で0の成分を削除した後、ALSを行う)
  
  ### 後処理
  # Cの各列をピッキングして、ピークがありそうな成分を選別する
  # matched filter関数を利用

  # ピーク領域と最大intensity領域が一致している場合はピークとみなす
  C_fin <- NULL;A_fin <- NULL
  for(i in 1:ncol(C)){
    if(which.max(matchedfilter(C[,i],0.5)[[1]])==which.max(C[,i]) && sum(A[i,]!=0)>=1){
      C_fin <- cbind(C_fin,C[,i])
      A_fin <- rbind(A_fin, A[i,])
    }
  }
  
  return(A_fin)
}
