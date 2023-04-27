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

### PCA
com <- 5
pca <- prcomp(Y)

### ICA
com <- 5
library(ica)
icacom <- ica(t(Y),com)

# ---------------------------------------------

### 関数化して、suppressする
C <- matrix(NA,com,nrow(Y))
for(k in 1:com){
  #S <- icacom$M[,k]
  S <- pca$x[,k]
  R <- cor(Y0,S)
  
  ## オリジナル
  if (R[which.max(abs(R))]<0){l
    S <- -S
  }
  S[S<0] <- 0
  C[k,] <- S
  
  # ## 両方の場合は両方にする(R>0.7)
  # if(R[which.max(abs(R))]>0.7){
  #   S1 <- S
  #   S1[S1<0] <- 0
  #   C <- rbind(C,S1)
  # }
  # if(R[which.min(abs(R))]<-0.7){
  #   S2 <- -S
  #   S2[S2<0] <- 0
  #   C <- rbind(C,S2)
  # }
}
# suppressWarnings(warning("testit"))

## NNLS
# library(nnls)
# A <- NULL
# for(i in 1:ncol(Y0)){
#   a <- nnls(t(C),Y0[,i])$x
#   A <- cbind(A,a)
# }

### マニュアル修正(2023.4.25)
# C[3,] <- -C[3,]+3000
# C[3,] <- C[3,]-mean(C[3,])
# index <- which(C[3,]<0)
# C[3,index] <- 0


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

# plot(mt/60,C[,2],type="l")
# par(new=T)
# plot(mt/60,C[,3],type="l")

### スペクトルマッチング
# スペクトルライブラリの読み込み
# マッチング

# Glutamine
# N-acetylcarnosine

library(MSnbase)
library(MsBackendMsp)

file_msp <- "C:/Users/yamamoto.HMT/Documents/R/msinfo/MSMS_Public_EXP_Pos_VS17.msp"
sp <- Spectra(file_msp, source = MsBackendMsp())

## 事前に一部のデータを書き出しておく
sp1 <- sp[19765]


metabolites <- sp@backend@spectraData@listData$name # 18237物質

#sum(metabolites=="Glutamine")
#sum(metabolites=="N-AcetylCarnosine")

### N-AcetylCarnosine
#19764	N-Acetylcarnosine; CE0; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19767	N-Acetylcarnosine; CE0; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19765	N-Acetylcarnosine; CE10; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19768	N-Acetylcarnosine; CE10; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19766	N-Acetylcarnosine; CE30; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19769	N-Acetylcarnosine; CE30; BKAYIFDRRZZKNF-VIFPVBQESA-N

index <- 19765 # CorrDecの結果(1番目の成分が最も相関が高い)
mz_target <- sp[index]@backend@spectraData@listData$mz[[1]]
int_target <- sp[index]@backend@spectraData@listData$intensity[[1]]

#plot(mz_target,int_target, type="h")

### Glutamine
# 1551	Glutamine
# 2404	Glutamine
# 4828	Glutamine
# 16573	GLUTAMINE
# 17146	GLUTAMINE
# 29971	Glutamine
# 29972	Glutamine
# 29973	Glutamine
# 29974	Glutamine
# 33380	Glutamine
# 84613	Glutamine (D)
# 84614	Glutamine (D)
# 84615	Glutamine (D)
# 84616	Glutamine (D)
# 84617	Glutamine (D)
# 84618	Glutamine (D)
# 84619	Glutamine (D)
# 84620	Glutamine (D)
# 84621	Glutamine (D)
# 84622	Glutamine (D)
# 84623	Glutamine (D)
# 84624	Glutamine (D)
# 90406	Glutamine (D)
# 84445	Glutamine (L)
# 84446	Glutamine (L)
# 84447	Glutamine (L)
# 84448	Glutamine (L)
# 84449	Glutamine (L)
# 84450	Glutamine (L)
# 84451	Glutamine (L)
# 84452	Glutamine (L)
# 84453	Glutamine (L)
# 84454	Glutamine (L)
# 84455	Glutamine (L)
# 84456	Glutamine (L)

index <- 29972 # 15eV (1番目の成分が最も相関が高い)
#index <- 29971 # 10eV

mz_target <- sp[index]@backend@spectraData@listData$mz[[1]]
int_target <- sp[index]@backend@spectraData@listData$intensity[[1]]

plot(mz_target,int_target, type="h")

# ------------------------------------------------------

# Metabolites
all <- cbind(mz_target,int_target)

#plot(all1,type="h", xlim=c(0,800))

# -------------------------------------

for(i in 1:com){
  
  mz_msp_i <- mz0
  int_msp_i <- A[i,]
  
  mz_msp_j <- all[,1]
  int_msp_j <- all[,2]
  
  s1 <- new("Spectrum2", mz=mz_msp_i, intensity=int_msp_i)
  s2 <- new("Spectrum2", mz=mz_msp_j, intensity=int_msp_j)
  
  r <- compareSpectra(s1, s2, fun="dotproduct") # 類似度尺度の検討
  print(r)
}
