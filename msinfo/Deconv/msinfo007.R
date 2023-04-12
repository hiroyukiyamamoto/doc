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

# EICを取得
y <- MSnbase::bin(swath_data_MS2_RT, binSize=0.01)

mz0 <- y[[1]]@mz

Y <- NULL;mt <- NULL
for(i in 1:length(swath_data_MS2_RT)){
  Y <- rbind(Y,y[[i]]@intensity)
  mt <- c(mt,y[[i]]@rt)
}

# 連続するデータポイントが3点以上
# ピークピッキング
index <- which(apply(Y,2,sd)!=0)
Z <- Y[,index]
mz <- mz0[index]

# W <- Z
# for(i in 1:ncol(Z)){
#   W[,i] <- Z[,i]/max(Z[,i])
# }

pca <- prcomp(Z)
#pca <- prcomp(Y[,index],scale=TRUE)

# R <- NULL
# for(i in 1:ncol(Z)){
#   R[i] <- cor.test(Z[,17],Z[,i])$estimate
# }
#
# index2 <- order(R, decreasing=TRUE)[1:10]
# -------------------------------------

plot(mt,pca$x[,1], type="l")
plot(mt,pca$x[,2], type="l")
plot(mt,pca$x[,3], type="l")




# -------------------------------------------
R <- cor(Z)
d <- as.dist(1 - R)
h <- hclust(d, method = "ave")
M <- cutree(h, k = 20)

# PC1のデータを生成
# 再度デコンボリューション

P <- NULL
for(i in 1:10){
  pca <- prcomp(Z[,which(M==i)])
  P <- cbind(P,pca$x[,1])
}

# M=3, M=7のplot
plot(P[,3], type="l")
plot(P[,7], type="l")


# a <- heatmap(R)

# 方法2
# スペクトル類似度を計算し、WGCNAでクラスタリングをして、
# クラスター中心からスペクトルを計算し、
# クロマトを逆行列で計算するのも良いかもしれない。

# 方法3
# 方法2と同じようにクロマトの類似度を計算して、クラスタリングする

# index <- c(17,13,11,1,5,4)
# plot(Z[,54], type="l")
#
# # -----------------------------------------
r <- NULL
for(i in 1:ncol(Z)){
 r[i] <- cor.test(Z[,135],Z[,i])$estimate
}
index1 <- which(r>0.8)

r1 <- NULL
for(i in 1:ncol(Z)){
 r1[i] <- cor.test(Z[,3],Z[,i])$estimate
}
index2 <- which(r1>0.8)

index <- c(index1,index2)

pca <- prcomp(Z[,index])




# ----------------------------------------

plot(mz[index],pca$rotation[,1], type="h", xlim=c(100,700))
plot(mz[index],pca$rotation[,2], type="h", xlim=c(100,700))
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
deconv <- osd(D=P, k=2, res.method='mcr')
deconv <- osd(D=Z, k=2, res.method='ica.osd')
deconv <- osd(D=Z, k=2, res.method='icr')

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
