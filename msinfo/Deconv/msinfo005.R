library(MSnbase)
library(loadings)

dia_file <- "C:/Users/hyama/Documents/msinfo/HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.mzML"
dia_data <- MSnbase::readMSData(dia_file, mode = "onDisk")

premz1 <- 300.1473
premz2 <- 290.1387
swath_data <- MSnbase::filterMsLevel(dia_data, msLevel=2L)
swath_data <- MSnbase::filterIsolationWindow(swath_data, mz=premz2)
swath_data_MS2_RT <- MSnbase::filterRt(swath_data, rt=c(170,190))
# MS1で時間幅を決める

# precursorに300.1473のみ設定している


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

plot(mz[index],pca$rotation[,1], type="l")
plot(mz[index],pca$rotation[,2], type="l")
plot(mz[index], pca$rotation[,3], type="l")

pc3 <- pca$rotation[,3]
pc3[pc3<0] <- 0

plot(mz[index], pc3, type="l")

plot(mz[index], pca$rotation[,2], type="l", xlim=c(100,200), ylim=c(0,1))

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
