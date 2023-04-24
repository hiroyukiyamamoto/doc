library(MSnbase)
library(loadings)

dia_file <- "C:/Users/hyama/Documents/msinfo/HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

premz <- 290.1387
swath_data <- MSnbase::filterMsLevel(dia_data, msLevel=2L)
swath_data <- MSnbase::filterIsolationWindow(swath_data, mz=premz)
swath_data_MS2_RT <- MSnbase::filterRt(swath_data, rt=c(170,190))
# MS1で時間幅を決める

y <- MSnbase::bin(swath_data_MS2_RT, binSize=0.01)

mz <- y[[1]]@mz

Y <- NULL;mt <- NULL
for(i in 1:length(swath_data_MS2_RT)){
  Y <- rbind(Y,y[[i]]@intensity)
  mt <- c(mt,y[[i]]@rt)
}
Y0 <- Y

### PCA

com <- 2
pca <- prcomp(Y)

# ---------------------------------------------

### 関数化して、suppressする
C <- NULL
for(k in 1:com){
  S <- pca$x[,k]
  R <- cor(Y0,S)

  ## 両方の場合は両方にする(R>0.7)
  if(R[which.max(abs(R))]>0.7){
    S1 <- S
    S1[S1<0] <- 0
    C <- rbind(C,S1)
  }
  if(R[which.min(abs(R))]<-0.7){
    S2 <- -S
    S2[S2<0] <- 0
    C <- rbind(C,S2)
  }
}
# suppressWarnings(warning("testit"))

## NNLS
# library(nnls)
# A <- NULL
# for(i in 1:ncol(Y0)){
#   a <- nnls(t(C),Y0[,i])$x
#   A <- cbind(A,a)
# }

# -------------------------------------

# number of component
lambda <- 0
maxiter <- 5

# initial value
C <- t(C)
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


plot(mt/60,C[,2],type="l")
par(new=T)
plot(mt/60,C[,3],type="l")






