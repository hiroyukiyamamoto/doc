# ---------
#   ICA
# ---------
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

# -------------------------------------

### ICA
library(ica)
com <- 5

icacom <- ica(t(Y0),com)

# ---------------------------------------------

### 関数化して、suppressする
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

### ALS
#library(osd)
#deconv <- osd(D=Y, k=com, res.method='mcr')

#plot(mt,deconv$C[,1])
#plot(deconv$C[,2])

#plot(mz0,deconv$S[,1],type="h")

# --------------------------------------------
#
# library(PSMatch)
#
# # GFSASSAR
# pepseq <- 'GFSASSAR'
# frag <- calculateFragments(pepseq, type=c("b", "y"))
# mz <- frag$mz
# int <- rep(1,length(mz))
#
# all1 <- cbind(mz,int)
#
# #plot(all1,type="h", xlim=c(0,800))
#
# # GFSANSAR
# pepseq <- 'GFSANSAR'
# frag <- calculateFragments(pepseq, type=c("b", "y"))
# mz <- frag$mz
# int <- rep(1,length(mz))
#
# all2 <- cbind(mz,int)
#
# #plot(all1,type="h", xlim=c(0,800))
#
# # library(wrProteo)
# # pep1 <- c(aa=pepseq)
# # convAASeq2mass(pep1, seqN=FALSE)
# #
# # convAASeq2mass(pep1, seqN=FALSE)/2 # 2価
#
# # -------------------------------------
# for(i in 1:com){
#
#   mz_msp_i <- mz0
#
#   S <- A[i,]
#
#   #int_msp_i <- deconv$S[,i]
#   int_msp_i <- S
#
#   mz_msp_j <- all1[,1]
#   int_msp_j <-all1[,2]
#
#   s1 <- new("Spectrum2", mz=mz_msp_i, intensity=int_msp_i)
#   s2 <- new("Spectrum2", mz=mz_msp_j, intensity=int_msp_j)
#
#   r <- compareSpectra(s1, s2, fun="dotproduct") # 類似度尺度の検討
#   print(r)
# }
#
# #plot(mz_msp_i,int_msp_i, type="h") # u_metabo[i]
# #plot(mz_msp_j,int_msp_j, type="h") # u_metabo[j]
# #
# # plot(mz_msp_i,int_msp_i/max(int_msp_i),type="h",xlim=c(min(mz_msp_i,mz_msp_j),1000),ylim=c(-1,1), lwd=2, xlab="m/z", ylab="")
# # par(new=T)
# # plot(mz_msp_j,-int_msp_j/max(int_msp_j),type="h",xlim=c(min(mz_msp_i,mz_msp_j),1000),ylim=c(-1,1), col="red",lwd=2, xlab="m/z", ylab="")
