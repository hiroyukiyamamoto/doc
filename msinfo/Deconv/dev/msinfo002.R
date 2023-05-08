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

