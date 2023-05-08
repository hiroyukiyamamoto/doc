library(MSnbase)
library(loadings)

dia_file <- "C:/Users/hyama/Documents/msinfo/HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

premz1 <- 300.1473 # Metochlopramide
premz2 <- 290.1387 # Norcocaine
swath_data <- MSnbase::filterMsLevel(dia_data, msLevel=2L)
swath_data <- MSnbase::filterIsolationWindow(swath_data, mz=premz2)
swath_data_MS2_RT <- MSnbase::filterRt(swath_data, rt=c(170,190))
# MS1で時間幅を決める

#10	SWATH	274	300

# precursorに300.1473のみ設定している

# Experiment	MS type	Min m/z	Max m/z
# 0	SCAN	50	1200
# 1	SWATH	50	75
# 2	SWATH	74	100
# 3	SWATH	99	125
# 4	SWATH	124	150
# 5	SWATH	149	175
# 6	SWATH	174	200
# 7	SWATH	199	225
# 8	SWATH	224	250
# 9	SWATH	249	275
# 10	SWATH	274	300
# 11	SWATH	299	325
# 12	SWATH	324	350
# 13	SWATH	349	375
# 14	SWATH	374	400
# 15	SWATH	399	425
# 16	SWATH	424	450
# 17	SWATH	449	475
# 18	SWATH	474	500


swath_data_MS1_RT <- MSnbase::filterRt(swath_data, rt=c(170,190))

# EICを取得

y <- MSnbase::bin(swath_data_MS2_RT, binSize=0.01)

mz <- y[[1]]@mz

Y <- NULL;mt <- NULL
for(i in 1:length(swath_data_MS2_RT)){
  Y <- rbind(Y,y[[i]]@intensity)
  mt <- c(mt,y[[i]]@rt)
}

library(ica)

icacom <- ica(t(Y), 5)

plot(icacom$M[,1], type="l")
plot(icacom$M[,2], type="l")
plot(icacom$M[,3], type="l")
plot(icacom$M[,4], type="l")
plot(icacom$M[,5], type="l")
plot(icacom$M[,6], type="l")
plot(icacom$M[,7], type="l")
plot(icacom$M[,8], type="l")
plot(icacom$M[,9], type="l")
plot(icacom$M[,10], type="l")
plot(icacom$M[,11], type="l")
plot(icacom$M[,12], type="l")
plot(icacom$M[,13], type="l")
plot(icacom$M[,14], type="l")
plot(icacom$M[,15], type="l")
plot(icacom$M[,16], type="l")
plot(icacom$M[,17], type="l")
plot(icacom$M[,18], type="l")
plot(icacom$M[,19], type="l")
plot(icacom$M[,20], type="l")



plot(icacom$S[1,], type="l")
plot(icacom$S[2,], type="l")
plot(icacom$S[3,], type="l")
plot(icacom$S[4,], type="l")
plot(icacom$S[5,], type="l")
plot(icacom$S[6,], type="l")
plot(icacom$S[7,], type="l")
plot(icacom$S[8,], type="l")
plot(icacom$S[9,], type="l")
plot(icacom$S[10,], type="l")


### pca

index <- which(apply(Y,2,sd)!=0)
Z <- Y[,index]

pca <- prcomp(Z)
#pca <- prcomp(Y[,index],scale=TRUE)

plot(mt/60,pca$x[,1], type="l")
plot(mt/60,pca$x[,2], type="l")
plot(pca$x[,3], type="l")

plot(mz[index],pca$rotation[,1], type="l")
plot(mz[index],pca$rotation[,2], type="l")
plot(mz[index], pca$rotation[,3], type="l")

pc3 <- pca$rotation[,3]
pc3[pc3<0] <- 0

plot(mz[index], pc3, type="l")

plot(mz[index], pca$rotation[,2], type="l", xlim=c(100,200), ylim=c(0,1))

summary(pca)

###
# 574

### Tsugawa
# 1216:184.015
# 1499:201.045
# 2018:227.055
# 574:136.075
# 1190:183.035
# 1767:212.035
# 276:105.035
# 1046:173.055
# 306:108.085

index <- c(1216,2018,574,1190,1767,276,1046,306)
plot(Z[,1216], type="l", ylim=c(0,10000))
for(i in 1:length(index)){
  par(new=T)
  plot(Z[,index[i]], type="l", ylim=c(0,10000))
}


plot(mt/60,Z[,index[3]], type="l", lty=2, ylim=c(0,3000)) # 136.075
par(new=T)
plot(mt/60,Z[,index[6]], type="l", lty=2, ylim=c(0,3000)) # 105.035
par(new=T)
plot(mt/60,Z[,index[8]], type="l", lty=2, ylim=c(0,3000)) # 108.085

par(new=T)
plot(mt/60,Z[,index[1]], type="l", ylim=c(0,3000)) # 136.075
par(new=T)
plot(mt/60,Z[,index[2]], type="l", ylim=c(0,3000)) # 136.075
par(new=T)
plot(mt/60,Z[,index[4]], type="l", ylim=c(0,3000)) # 136.075
par(new=T)
plot(mt/60,Z[,index[5]], type="l", ylim=c(0,3000)) # 136.075


#
# library(osd)
# deconv <- osd(D=Y[,index], k=5, res.method='icr')
# # icrとica.osdは時間が掛るのでmcrを使う
#
# par(mfrow = c(2, 1))
# plot(mt,deconv$C[,1],type="l")
# plot(mt,deconv$C[,2],type="l")
# plot(mt,deconv$C[,3],type="l")
# plot(mt,deconv$C[,4],type="l")
# plot(mt,deconv$C[,5],type="l")

# plot(mz[index],deconv$S[,1], type="h")
# plot(mz[index],deconv$S[,2], type="h")
