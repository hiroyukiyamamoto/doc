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
# com <- 5
# pca <- prcomp(Y)

### ICA
com <- 5
library(ica)
icacom <- ica(t(Y),com)

A <- icacom$S

# ----------------------------------------------

library(MSnbase)
library(MsBackendMsp)

#file_msp <- "C:/Users/yamamoto.HMT/Documents/R/msinfo/MSMS_Public_EXP_Pos_VS17.msp"
#sp <- Spectra(file_msp, source = MsBackendMsp())
# 
# metabolites <- sp@backend@spectraData@listData$name # 18237物質

#sum(metabolites=="Glutamine")
#sum(metabolites=="N-AcetylCarnosine")

### N-AcetylCarnosine
#19764	N-Acetylcarnosine; CE0; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19767	N-Acetylcarnosine; CE0; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19765	N-Acetylcarnosine; CE10; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19768	N-Acetylcarnosine; CE10; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19766	N-Acetylcarnosine; CE30; BKAYIFDRRZZKNF-VIFPVBQESA-N
#19769	N-Acetylcarnosine; CE30; BKAYIFDRRZZKNF-VIFPVBQESA-N

#index <- 19764
#index <- 19765 # CorrDecの結果
#index <- 19766 # 気になる
#index <- 19767
#index <- 19768 # 気になる
#index <- 19769
#mz_target <- sp[index]@backend@spectraData@listData$mz[[1]]
#int_target <- sp[index]@backend@spectraData@listData$intensity[[1]]

# --------------------------------
#   N-AcetylCarnosine (CorrDec)
# --------------------------------
#NAC <- data.frame(
#  mz_target = sp[index]@backend@spectraData@listData$mz[[1]],
#  int_target = sp[index]@backend@spectraData@listData$intensity[[1]]
#)

# save(NAC,"C:/Users/yamamoto.HMT/Documents/R/msinfo/CorrDec/NAC.RData")

load("C:/Users/yamamoto.HMT/Documents/R/msinfo/CorrDec/NAC.RData")

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

#index <- 29972 # 15eV (1番目の成分)

# index <- 29972 # 10eV
# mz_target <- sp[index]@backend@spectraData@listData$mz[[1]]
# int_target <- sp[index]@backend@spectraData@listData$intensity[[1]]
# 
# plot(mz_target,int_target, type="h")

# ------------------------------------------------------

# Metabolites
all <- cbind(NAC$mz_target,NAC$int_target)

#plot(all1,type="h", xlim=c(0,800))

# -------------------------------------

for(i in 1:com){
  
  A1 <- A[,i]
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






