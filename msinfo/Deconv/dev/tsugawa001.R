rm(list=ls(all=TRUE))

### MS1でピッキングして、ピーク領域の範囲を決めた後で、MS2のデコンボリューションを行う
### 成分数を5としてデコンボリューションして、ALS5回繰り返したもので良しとする
### その後、スペクトルライブラリを使って類似度を計算し、アノテーションをあてる
### MS1のアノテーションがMS2を加えることによって改善するかを確認する

source("C:/R/function/errppm.R")
source("C:/Users/yamamoto.HMT/Documents/R/msinfo/DIA/deconvICA.R")

library(MSnbase)
library(loadings)
library(xcms)

dia_file <- "C:/Users/yamamoto.HMT/Documents/R/msinfo/HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.mzML"
swath_data <- readMSData(dia_file, mode = "onDisk")
swath_data0 <- swath_data

table(msLevel(swath_data)) # 1234, 22212

### isolation windowの設定確認 ----------------------------
head(fData(swath_data)[, c("isolationWindowTargetMZ",
                           "isolationWindowLowerOffset",
                           "isolationWindowUpperOffset",
                           "msLevel", "retentionTime")])

head(isolationWindowUpperMz(swath_data))
head(isolationWindowLowerMz(swath_data))

table(isolationWindowTargetMz(swath_data))

# ---------------------------------------------------------

### MS/MSスペクトルライブラリ読み込み ---------------------
library(MsBackendMsp)
file_msp <- "C:/Users/yamamoto.HMT/Documents/R/msinfo/MSMS_Public_EXP_Pos_VS17.msp" # 72439 spectra
sp <- Spectra(file_msp, source = MsBackendMsp())
sp_premz <- sp@backend@spectraData$precursorMz # precursorのm/z

# ---------------------------------------------------------

### MS1 ピッキング
cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3, 30))
swath_data <- findChromPeaks(swath_data, param = cwp)

### ピークリスト
peaklist <- chromPeaks(swath_data) # 1327ピーク

### 
N <- NULL; M <- NULL; R <- NULL; Q <- NULL
for(i in 1:nrow(peaklist)){
  print(i)
  
  targetmz <- peaklist[i,1]
  targetRT <- peaklist[i,4]
  
  # isolation windowで選別する
  swath_data_iwindow <- filterIsolationWindow(swath_data, mz = targetmz)

  rtmin <- peaklist[i,5]
  rtmax <- peaklist[i,6]
  
  # MS1のピーク領域(RT)とそのMS2
  #index <- which((rtmin < targetRT) & (targetRT < rtmax))
  #spec_all <- swath_data_iwindow[index]
  
  # MS2のmatrixに整形
  swath_data_MS2_RT <- MSnbase::filterRt(swath_data_iwindow, rt=c(rtmin,rtmax))
  
  ## MS/MSのサイズ
  N[i] <- length(swath_data_MS2_RT)
  
  # ICAの成分数が5なので、最低でも5個のMS/MSスペクトルが必要
  if(length(swath_data_MS2_RT)>=5){
    y <- MSnbase::bin(swath_data_MS2_RT, binSize=0.01)
    
    mz <- y[[1]]@mz
    Y <- NULL;mt <- NULL
    for(j in 1:length(swath_data_MS2_RT)){
      Y <- rbind(Y,y[[j]]@intensity)
      mt <- c(mt,y[[j]]@rt)
    }
    Y0 <- Y
    
    # デコンボリューション
    A <- deconvICA(Y0)
  }
  
  # デコンボリューションしないときは、RTが最も近い時のMS/MSスペクトルを選ぶ
  if(length(swath_data_MS2_RT)<5){
    if (N[i] >= 1){
      index_min <- which.min(abs(swath_data_MS2_RT@featureData@data$retentionTime-targetRT))
      MS2_sample <- swath_data_MS2_RT[[index_min]]
    }
  }
  
  ## -----------------------------
  ##　デコンボリューションあり
  ## -----------------------------
  if(length(swath_data_MS2_RT)>=5){
    
    # MS/MSスペクトルマッチング
    index <- which(abs(e_ppm(targetmz,sp_premz))<5)
  
    # 類似度計算
    q <- NULL
    for(com in 1:nrow(A)){
      spec <- A[com,]
      r <- NULL
      if(length(index)>=1){
        for(k in 1:length(index)){
          # スペクトルライブラリ
          mz_msp <- sp[index[k]]@backend@spectraData@listData$mz[[1]]
          int_msp <- sp[index[k]]@backend@spectraData@listData$intensity[[1]]
        
          # 実測スペクトル
          mz_sample <- mz[which(spec!=0)]
          int_sample <- spec[which(spec!=0)]
        
          s1 <- new("Spectrum2", mz=mz_msp, intensity=int_msp)
          s2 <- new("Spectrum2", mz=mz_sample, intensity=int_sample)
        
          r[k] <- compareSpectra(s1, s2, fun="dotproduct") # 類似度尺度の検討
        }
      }
      q[com] <- max(r)
    }
    Q[i] <- max(q)
  }else{

    ## ----------------------------
    ##　デコンボリューションなし
    ## ----------------------------
    # 類似度の分布を比較する(デコンボリューションありなし)
  
    # MS/MSスペクトルマッチング
    index <- which(abs(e_ppm(targetmz,sp_premz))<5)

    # 類似度計算
    ##  類似度は、デコンボリューションした時としていない時の最大のものを採用すれば、
    ## トラブルは起きにくいか？確認する
    r <- NULL
    if(length(index)>=1){
      for(k in 1:length(index)){
        # スペクトルライブラリ
        mz_msp <- sp[index[k]]@backend@spectraData@listData$mz[[1]]
        int_msp <- sp[index[k]]@backend@spectraData@listData$intensity[[1]]
      
        # 実測スペクトル
        mz_sample <- MS2_sample@mz
        int_sample <- MS2_sample@intensity

        s1 <- new("Spectrum2", mz=mz_msp, intensity=int_msp)
        s2 <- new("Spectrum2", mz=mz_sample, intensity=int_sample)
      
        r[k] <- compareSpectra(s1, s2, fun="dotproduct") # 類似度尺度の検討
      }
      R[i] <- max(r)
    }else{
      R[i] <- NA
    }
  }
}

hist(R[!is.na(R)], breaks=100)

Q[Q==-Inf] <- NA
hist(Q[!is.na(Q)], breaks=100)

