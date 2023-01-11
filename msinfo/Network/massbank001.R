rm(list=ls(all=TRUE))

# library(mzR)
library(MsBackendMsp)
library(MSnbase)

#file_msp <- "C:/Users/yamamoto/Documents/R/msinfo/Network/MoNA-export-MassBank.msp" # 72439 spectra
#sp <- Spectra(file_msp, source = MsBackendMsp())

#save(sp, file="C:/Users/yamamoto/Documents/R/msinfo/Network/sp_MassBank.RData")
load(file="C:/Users/yamamoto/Documents/R/msinfo/Network/sp_MassBank.RData")

# 物質名一覧
metabolites <- sp@backend@spectraData@listData$Name # 18237物質

# 同じmetabolitesのMS/MSスペクトルを統合する
u_metabo <- unique(metabolites)
# 
# spall <- NULL
# for(i in 1:length(u_metabo)){
#   print(i)
#   index <- which(u_metabo[i]==metabolites)
#   if(length(index)>=2){
#     sp1 <- NULL
#     for(j in 1:length(index)){
#       mz <- sp@backend@spectraData@listData$mz[[index[j]]]
#       intensity <- sp@backend@spectraData@listData$intensity[[index[j]]]
#       sp0 <- new("Spectrum1", mz = mz, intensity = intensity)
#       sp1 <- c(sp1,sp0)
#     }
#     spctra <- MSpectra(sp1)
#     csp <- combineSpectra(spctra, mzd = 0.05, intensityFun = max)@listData$`1` # 統合するパラメーターは要検討
#     spall[[i]] <- csp
#   }else{
#     mz <- sp@backend@spectraData@listData$mz[[index]]
#     intensity <- sp@backend@spectraData@listData$intensity[[index]]
#     sp0 <- new("Spectrum1", mz = mz, intensity = intensity)
#     spall[[i]] <- sp0
#   }
# }

# save(spall, file="C:/Users/yamamoto/Documents/R/msinfo/Network/spall_metabolite.RData")
  
# --------------------------------------------------------------

load(file="C:/Users/yamamoto/Documents/R/msinfo/Network/spall_metabolite.RData")

# # 類似度計算
# R <- matrix(0,length(u_metabo),length(u_metabo))
# for(i in 1:length(u_metabo)){
#   print(i)
#   r <- NULL
#   for(j in 1:length(u_metabo)){
#     mz_msp_i <- spall[[i]]@mz
#     int_msp_i <- spall[[i]]@intensity
#     mz_msp_j <- spall[[j]]@mz
#     int_msp_j <- spall[[j]]@intensity
#     
#     s1 <- new("Spectrum2", mz=mz_msp_i, intensity=int_msp_i)
#     s2 <- new("Spectrum2", mz=mz_msp_j, intensity=int_msp_j)
#     
#     r[j] <- compareSpectra(s1, s2, fun="dotproduct") # 類似度尺度の検討
#   }
#   R[i,] <- r
# }

# ------------------------------------------------------

# 100物質 ランダムサンプリング
#index <- sample(c(1:length(u_metabo)))[1:100]
#save(index, file="C:/Users/yamamoto/Documents/R/msinfo/Network/100.RData")
load(file="C:/Users/yamamoto/Documents/R/msinfo/Network/100.RData")

# 類似度計算
R <- NULL
for(i in index){
  print(i)
  r <- NULL;k <- 1
  for(j in index){
    mz_msp_i <- spall[[i]]@mz
    int_msp_i <- spall[[i]]@intensity
    mz_msp_j <- spall[[j]]@mz
    int_msp_j <- spall[[j]]@intensity
    
    s1 <- new("Spectrum2", mz=mz_msp_i, intensity=int_msp_i)
    s2 <- new("Spectrum2", mz=mz_msp_j, intensity=int_msp_j)
    
    r[k] <- compareSpectra(s1, s2, fun="dotproduct") # 類似度尺度の検討
    k <- k+1
  }
  R <- rbind(R,r)
}

#row.names(R) <- index
#colnames(R) <- index
#
#heatmap(R,Colv = NA, Rowv = NA)

th <- 0.5
A <- abs(R)
A[R >= th] <- 1
A[R < th] <- 0

g <- graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected", diag = F )

### ネットワークを生成
plot(g, layout = layout_with_mds, vertex.label.cex = 1, vertex.frame.color = "white", vertex.size = 5, vertex.color = "lightblue")

# ------------------------------------------------------

# 6991と7560のスペクトルの類似度が最も高い
# 実際に表示してみる

i <- 6991
mz_msp_i <- spall[[i]]@mz
int_msp_i <- spall[[i]]@intensity

j <- 7560
mz_msp_j <- spall[[j]]@mz
int_msp_j <- spall[[j]]@intensity

plot(mz_msp_i,int_msp_i, type="h") # u_metabo[i]
plot(mz_msp_j,int_msp_j, type="h") # u_metabo[j]

plot(mz_msp_i,int_msp_i/max(int_msp_i),type="h",xlim=c(min(mz_msp_i,mz_msp_j),max(int_msp_i,int_msp_j)),ylim=c(-1,1), lwd=2, xlab="m/z", ylab="")
par(new=T)
plot(mz_msp_j,-int_msp_j/max(int_msp_j),type="h",xlim=c(min(mz_msp_i,mz_msp_j),max(int_msp_i,int_msp_j)),ylim=c(-1,1), col="red",lwd=2, xlab="m/z", ylab="")
