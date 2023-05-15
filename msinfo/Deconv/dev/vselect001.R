rm(list=ls(all=TRUE))

# -----------------
#   CorrDec (PCA)
# -----------------
library(MSnbase)
library(loadings)

dia_file <- "C:/Users/yamamoto.HMT/Documents/R/msinfo/CorrDec/QC1.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

swath_data <- MSnbase::filterMsLevel(dia_data, msLevel=2L)
swath_data_MS2_RT <- MSnbase::filterRt(swath_data, rt=c(8.8*60,9.1*60)) # 8.8～9.1min

y <- MSnbase::bin(swath_data_MS2_RT, binSize=0.01)

mz <- y[[1]]@mz
mz0 <- mz

Y <- NULL;mt <- NULL
for(i in 1:length(swath_data_MS2_RT)){
  if(swath_data_MS2_RT[[i]]@collisionEnergy==10){
    Y <- rbind(Y,y[[i]]@intensity)
    mt <- c(mt,y[[i]]@rt)
  }
}
Y0 <- Y

# -------------------
#   変数選択(手動)
# -------------------
### binSize=0.01の時
# 147.077: 10708
# 130.050: 9005 (クロマトが無い) 9006にある
# 84.045: 4405
# 56.049: 1605

# 269.125: 22913
# 223.119: 18312
# 156.077: 11608
# 110.071: 7008

### binSize=0.001の時
# 147.077: 107077
# 130.050: 90050 (クロマトが無い) 90051にある
# 84.045: 44045
# 56.049: 16049

# 269.125: 229125
# 223.119: 183119
# 156.077: 116077
# 110.071: 70071

which.min(abs(mz-147.077))

### binsize=0.01
index <- c(10708,9006,4405,1605,22913,18312,11608,7008)
### binsize=0.001
#index <- c(107077,90051,44045,16049,229125,183119,116077,70071)

Z <- Y0[,index]

#Z <- Y0[,apply(Y,2,sd)!=0] # 変数削減
### MS2クロマトをピークピッキングして変数を減らすのもありかもしれない

# ---------------------------------------------

### ICA
com <- 5
library(ica)
icacom <- ica(t(Z),com)
# icacom <- ica(Z,com)

A <- icacom$S

plot(icacom$M[,1], type="l")
plot(icacom$S[,1], type="h")

plot(icacom$M[,2], type="l")
plot(icacom$S[,2], type="h")

plot(icacom$M[,3], type="l")
plot(icacom$S[,3], type="h")

plot(icacom$M[,4], type="l")
plot(icacom$S[,4], type="h")

plot(icacom$M[,5], type="l")
plot(icacom$S[,5], type="h")

# -------------------------------------

### 関数化して、suppressする
com <- 5
C <- matrix(NA,com,nrow(Z))
for(k in 1:com){
  S <- icacom$M[,k]
  R <- cor(Z,S)
  
  if (R[which.max(abs(R))]<0){
    S <- -S
  }
  S[S<0] <- 0
  C[k,] <- S
}

# -------------------------------------

# number of component
lambda <- 0
maxiter <- 5

C <- t(C)
X <- Z

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
