# ----------------------
#   ガウスフィルタ関数
# ----------------------
matchedfilter <- function (xt,fwhm,method="gaussian"){
  
  sigma <- fwhm/2.3548
  
  # 初期値
  M <- length(xt)
  Ybis <- as.vector(matrix(0,1,M))    # 全部0(端)
  
  # 設定値
  m <- ceiling(4*sigma)
  
  # ガウシアンフィルタ
  w　<- -m:m
  fbis <- as.vector((1-w^2/sigma^2)*exp(-w^2/2/sigma^2))    # ガウシアンフィルタ
  
  # 2次微分
  for (tx in (m+1):(M-(m+1)) ){
    if (tx-m >= 1 & tx+m <= M){    
      Ybis[tx] <- sum(xt[c((tx-m):(tx+m))]*(fbis))/(sqrt(sum(fbis*fbis)))
    }
  }
  Y <- list(Ybis,fbis)
  return(Y)
}

# -------------
#   XCMSのfft
# -------------
filtfft <- function(y, filt) {
  
  yfilt <- numeric(length(filt))
  yfilt[1:length(y)] <- y
  yfilt <- fft(fft(yfilt, inverse = TRUE) * filt)
  
  Re(yfilt[1:length(y)])
}
