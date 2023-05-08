library(xcms)

dia_file <- "C:/Users/hyama/Documents/msinfo/HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.mzML"
dia_data <- readMSData(dia_file, mode = "onDisk")

x1 <- xcms::chromatogram(dia_data, mz=c(300.1,300.2),msLevel = 1L, aggregationFun = "sum")
x2 <- xcms::chromatogram(dia_data, mz=c(290.1,290.2),msLevel = 1L, aggregationFun = "sum")

plot(x1@.Data[[1]]@rtime/60,x1@.Data[[1]]@intensity, type="l", xlim=c(2.8,3.2), ylim=c(0,10000))
par(new=T)
plot(x2@.Data[[1]]@rtime/60,x2@.Data[[1]]@intensity, type="l", col="red", xlim=c(2.8,3.2), ylim=c(0,10000))

plot(x1@.Data[[1]]@rtime,x1@.Data[[1]]@intensity, type="l", xlim=c(172,185), ylim=c(0,10000))
par(new=T)
plot(x2@.Data[[1]]@rtime,x2@.Data[[1]]@intensity, type="l", col="red", xlim=c(172,185), ylim=c(0,10000))


library(mzR)

aa <- openMSfile(dia_file)

k <- runInfo(aa)$scanCount

all <- NULL
for(i in 1:k){
  all[[i]] <- peaks(aa,i)
  header(aa,i)
}

# mzRでデータを読み込んで整理するか、
# MSnBaseで読み込んで整理するのが良いか、確認する
