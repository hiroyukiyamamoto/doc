library(MSnbase)
library(loadings)

dia_file <- "C:/Users/hyama/Documents/msinfo/CN20161108_SAM_SPECTER_NB2p81_11.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

premz1 <- 391.6932
premz2 <- 405.1987
swath_data <- MSnbase::filterMsLevel(dia_data, msLevel=2L)
swath_data <- MSnbase::filterIsolationWindow(swath_data, mz=premz1)
swath_data_MS2_RT <- MSnbase::filterRt(swath_data, rt=c(0,500))

# isolation windowをどれを使っているか確認する必要がある
# 400.4369 422.4369
# 389.4319 411.4319

# premz1を使うと、デコンボリューションがうまくいかない
# おそらくMSMSスペクトルがデコンボリューションに合わないデータになっている
# premz2を使って400.4369 422.4369

# MS1で時間幅を決める

#swath_data_MS1_RT <- MSnbase::filterRt(swath_data1, rt=c(170,190))

# EICを取得

y <- MSnbase::bin(swath_data_MS2_RT, binSize=0.01)

mz <- y[[1]]@mz

Y <- NULL;mt <- NULL
for(i in 1:length(swath_data_MS2_RT)){
  Y <- rbind(Y,y[[i]]@intensity)
  mt <- c(mt,y[[i]]@rt)
}

# Y <- Y[seq(1,nrow(Y),2),]

### pca

 N <- NULL
 for(i in 1:ncol(Y)){
   y <- Y[,i]
   index0 <- which(y>0)
   # 連続して5ポイント存在するピーク
   N[i] <- sum(diff(index0) %in% c(1,1,1,1,1,1,1,1,1,1))
 }

# index <- which(N>=10)
# 連続するデータポイントが3点以上
# ピークピッキング
index <- which(apply(Y,2,sd)!=0)
Z <- Y[,index]

pca <- prcomp(Z)
#pca <- prcomp(Y[,index],scale=TRUE)

plot(Z[,1], type="l")

# R <- NULL
# for(i in 1:ncol(Z)){
#   R[i] <- cor.test(Z[,17],Z[,i])$estimate
# }
#
# index2 <- order(R, decreasing=TRUE)[1:10]

ospca2 <- function(X, D, kappa = 0.999, M = diag(1, nrow(X))){
  MX <- scale(M %*% X, scale=FALSE)
  X <- scale(X, scale=FALSE)
  E <- (1 - kappa) * diag(1, ncol(MX)) + kappa * t(MX) %*%
    t(D) %*% D %*% MX
  G <- chol(solve(E))  # かなり時間が掛かって計算できない。
  # カーネルにすると早くなるはず
  # 片側カーネル平滑化を考えるのが良い
  # 計算できたとしても、上手くいく保証はない
  # 変数をかなり削れば計算できるはず→削ったが計算できなかった

  # ただし、平滑化os-pcaを高速化すれば解決する問題かどうかが不明

  W0 <- svd(G %*% t(MX) %*% MX)$v
  t <- X %*% W0
  Mt <- MX %*% W0
  R <- chol(E)
  z <- svd(t(MX) %*% MX %*% solve(R))$v
  W2 <- solve(R) %*% z
  Ms <- MX %*% W2
  list(P = W0, T = t, MT = Mt, Q = W2, U = Ms)
}

class <- c(1:103)
D0 <- factor(class)
D1 <- model.matrix(~ D0 + 0)
D <- diff(D1)

ospca <- ospca2(Z,D, kappa=0.0000001)
# OS-PCAの高速化
# -------------------------------------

# plot(ospca$T[,1], type="l")
# par(new=T)
# plot(ospca$T[,2], type="l")
# plot(ospca$T[,3], type="l")



# -------------------------------------

plot(mt,pca$x[,1], type="l")
plot(mt,pca$x[,2], type="l")
plot(mt,pca$x[,3], type="l")


plot(mt,pca$x[,4], type="l")
plot(mt,pca$x[,5], type="l")
plot(mt,pca$x[,6], type="l")
plot(mt,pca$x[,7], type="l")
plot(mt,pca$x[,8], type="l")
plot(pca$x[,9], type="l")
plot(pca$x[,10], type="l")
plot(pca$x[,11], type="l")
plot(pca$x[,12], type="l")
plot(pca$x[,13], type="l")
plot(pca$x[,14], type="l")
plot(pca$x[,15], type="l")
plot(pca$x[,16], type="l")
plot(pca$x[,17], type="l")
plot(pca$x[,18], type="l")
plot(pca$x[,19], type="l")
plot(pca$x[,20], type="l")

plot(mz[index],pca$rotation[,1], type="l", xlim=c(100,700))
plot(mz[index],pca$rotation[,2], type="l", xlim=c(100,700))
plot(mz[index], pca$rotation[,3], type="l", xlim=c(100,700))
plot(mz[index], pca$rotation[,4], type="l", xlim=c(100,700))
plot(mz[index], pca$rotation[,5], type="l", xlim=c(100,700))
plot(mz[index], pca$rotation[,6], type="l", xlim=c(100,700))

pc3 <- pca$rotation[,3]
pc3[pc3<0] <- 0

plot(mz[index], pc3, type="l")

plot(mz[index], pca$rotation[,2], type="l", xlim=c(100,200), ylim=c(0,1))

summary(pca)

# --------------------------------------------

#
library(osd)
deconv <- osd(D=Z, k=5, res.method='mcr')
# # icrとica.osdは時間が掛るのでmcrを使う
#
# par(mfrow = c(2, 1))
plot(mt,deconv$C[,1],type="l")
plot(mt,deconv$C[,2],type="l")


# plot(mz[index],deconv$S[,1], type="h")
# plot(mz[index],deconv$S[,2], type="h")

# --------------------------------------------

library(PSMatch)

# GFSASSAR
pepseq <- 'GFSASSAR'
frag <- calculateFragments(pepseq, type=c("b", "y"))
mz <- frag$mz
int <- rep(1,length(mz))

all1 <- cbind(mz,int)

plot(all1,type="h", xlim=c(0,800))


library(wrProteo)
pep1 <- c(aa=pepseq)
convAASeq2mass(pep1, seqN=FALSE)

convAASeq2mass(pep1, seqN=FALSE)/2 # 2価
