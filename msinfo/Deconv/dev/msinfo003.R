library(MSnbase)

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

###
library(osd)
deconv <- osd(D=Y, k=5, res.method='mcr')
# icrとica.osdは時間が掛るのでmcrを使う


par(mfrow = c(2, 1))
plot(mt,deconv$C[,1],type="l")
plot(mt,deconv$C[,2],type="l")



#par(mfrow = c(3, 2))
#plotOSDres(deconv, type='eic',1)
#plotOSDres(deconv, type='s',1)
#plotOSDres(deconv, type='eic',2)
#plotOSDres(deconv, type='s',2)

