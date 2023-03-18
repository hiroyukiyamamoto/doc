library(MSnbase)
library(loadings)

dia_file <- "C:/Users/hyama/Documents/msinfo/HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

premz <- 300.1473
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

### pca

index <- which(apply(Y,2,sd)!=0)

pca <- prcomp(Y[,index])
#pca <- prcomp(Y[,index],scale=TRUE)

plot(pca$x[,1], type="l")
plot(pca$x[,2], type="l")
plot(pca$x[,3], type="l")

pca <- prcomp(t(Y[,index]))
#pca <- prcomp(t(Y[,index]),scale=TRUE)

plot(pca$rotation[,1], type="l")
par(new=T)
plot(-pca$rotation[,2], type="l")
plot(pca$rotation[,3], type="l")

### pls

# class <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10)
class <- c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6)

Y0 <- factor(class)
Y1 <- model.matrix(~ Y0 + 0)

library(chemometrics)
pls <- pls_eigen(scale(Y[,index],scale=FALSE),scale(Y1,scale=FALSE),5)

plot(pls$T[,1], type="l")
plot(pls$T[,2], type="l")
plot(pls$T[,3], type="l")

pls <- pls_eigen(scale(Y[,index]),scale(Y1),5)

plot(pls$T[,1], type="l")
plot(pls$T[,2], type="l")
plot(pls$T[,3], type="l")

pls <- pls_eigen(Y[,index],Y1,5)

plot(-pls$T[,1], type="l")
par(new=T)
plot(pls$T[,2], type="l")
plot(pls$T[,3])

### os-pca

ospca2 <- function(X, D, kappa = 0.999, M = diag(1, nrow(X))){
  MX <- scale(M %*% X, scale=FALSE)
  X <- scale(X, scale=FALSE)
  E <- (1 - kappa) * diag(1, ncol(MX)) + kappa * t(MX) %*%
    t(D) %*% D %*% MX
  G <- chol(solve(E))
  W0 <- svd(G %*% t(MX) %*% MX)$v
  t <- X %*% W0
  Mt <- MX %*% W0
  R <- chol(E)
  z <- svd(t(MX) %*% MX %*% solve(R))$v
  W2 <- solve(R) %*% z
  Ms <- MX %*% W2
  list(P = W0, T = t, MT = Mt, Q = W2, U = Ms)
}

class <- c(1:29)
D0 <- factor(class)
D1 <- model.matrix(~ D0 + 0)
D <- diff(D1)

ospca <- ospca2(Y[,index],D, kappa=0.001)

plot(ospca$T[,1], type="l")
par(new=T)
plot(ospca$T[,2], type="l")
plot(ospca$T[,3], type="l")

plot(ospca$U[,1], type="l")
par(new=T)
plot(ospca$U[,2], type="l")
plot(ospca$U[,3], type="l")

# Uを考えると結構大変なので、とりあえずは考えない
