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

# --------------------
#   Non-negative ICA
# --------------------
# 完全にnon-negativeにはならないか？
# 実装が間違っていないか、理論を再度確認する必要がある
### スペクトル同士を独立として計算する場合になっているか？

com <- 5
W <- diag(1,com)

### オリジナルのプログラム
#usv <- svd(Y)
#M <- t(usv$v[,1:com])

### PCA
pca <- prcomp(Y)
M <- t(pca$rotation[,c(1:com)])


X <- M

library(expm)

eta <- 0.01
for(h in 1:1000){
  Y <- X
  
  Yp <- Y
  Yp[Yp < 0] <- 0
  
  Yn <- Y
  Yn[Yn > 0] <- 0
  
  W <- W+expm(-eta*(Yn%*%t(Yp)-Yp%*%t(Yn)))
  
  X <- W%*%M
  
}

X[X<0] <- 0 # positive

### 
# X: 2×42400
# Z : 29×42400

#C <- t(solve(X%*%t(X))%*%X%*%t(Z)) # 29×2
Y <- Y0
# C <- t(solve(X%*%t(X))%*%X%*%t(Y)) # 29×2

library(nnls)

C <- NULL
for(i in 1:nrow(Y)){
  c0 <- nnls(t(X),(Y[i,]))$x
  C <- cbind(C,c0)
}
C <- t(C)

CnnICA <- C

# ----------------
#   ALS
# ----------------

# number of component
lambda <- 0.01
maxiter <- 10

# initial value
X <- Y0

E <- NULL
for (k in 1:maxiter){
  A <- solve(t(C)%*%C+lambda*diag(1,ncol(C)))%*%t(C)%*%X
  A[A<0] <- 0
  
  # normalized constraint
  #for (i in 1:com){
  #  A[i,] <- A[i,]/as.numeric(sqrt(t(A[i,])%*%A[i,]))
  #}
  
  C <- X%*%t(A)%*%solve(A%*%t(A)+lambda*diag(1,ncol(C)))
  C[C<0] <- 0
  
  # error
  E[k] <- norm(X-C%*%A);
}


# -------------------------
#   スペクトルマッチング
# -------------------------

load("C:/Users/yamamoto.HMT/Documents/R/msinfo/CorrDec/NAC.RData")

# Metabolites
all <- cbind(NAC$mz_target,NAC$int_target)

#plot(all1,type="h", xlim=c(0,800))

for(i in 1:com){
  
  A1 <- A[i,]
  A1[A1 < 0] <- 0
  
  mz_msp_i <- mz0
  int_msp_i <- A1
  
  mz_msp_j <- all[,1]
  int_msp_j <- all[,2]
  
  s1 <- new("Spectrum2", mz=mz_msp_i, intensity=int_msp_i)
  s2 <- new("Spectrum2", mz=mz_msp_j, intensity=int_msp_j)
  
  r <- compareSpectra(s1, s2, fun="dotproduct") # 類似度尺度の検討
  print(r)
}




