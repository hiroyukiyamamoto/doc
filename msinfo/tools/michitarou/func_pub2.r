# -----------------
# 	PubChem検索
# -----------------
func_pub2 <- function(mf,preID){

  library(rcdk)
  
  # ---------------------
  #   ファイル一覧削除
  # ---------------------
  foldername <- "C:/R/sdf/"
  flist <- paste0(foldername,list.files(foldername))
  
  if (length(list.files(foldername))>0){
    tt <- file.remove(flist)  
  }
  
  # -----------------------
  #   分子式から構造検索
  # -----------------------
  # mf <- "C7H12N2O5" # 未知ピークの分子式
  
  # PubChemで分子式検索
  url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/",mf,"/cids/TXT",sep="")
  
  # CIDのリストをダウンロード
  download.file(url, "C:/R/pub1.txt", method = "auto", quiet=FALSE)
  
  # -------------------------------
  #   sdfファイルのダウンロード
  # -------------------------------
  # PubChemで分子式検索
  x <- read.csv("C:/R/pub1.txt", header=FALSE)
  CID <- x[[1]]
  
  FP <- NULL
  for(i in 1:nrow(x)){
    url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",CID[i],"/SDF",sep="")
    download.file(url, paste0("C:/R/sdf/", CID[i], ".sdf"), method = "auto", quiet=FALSE)
    mols <- load.molecules(paste0("C:/R/sdf/",CID[i],".sdf"))
    fp <- get.fingerprint(mols[[1]], type="pubchem")
    
    z <- as.numeric(matrix(0,881)) # 初期値
    z[fp@bits] <- 1
    
    FP <- rbind(FP,z)
  }
  
  # -------------------------------
  #   前駆体のpubchem fingerprint
  # -------------------------------
  # pubchem CID（化合物リスト）
  # preID <- 119 # 要設定
  
  # ダウンロード
  url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",preID,"/SDF",sep="")
  download.file(url, "C:/R/pub.sdf", method = "auto", quiet=FALSE)
  mols <- load.molecules('C:/R/pub.sdf')
  fp <- get.fingerprint(mols[[1]], type="pubchem")
  
  preFP <- as.numeric(matrix(0,881)) # 初期値
  preFP[fp@bits] <- 1
  
  # -----------------
  #   Tanimoto係数
  # -----------------
  R <- NULL
  for(i in 1:nrow(FP)){
    
    a11 <- sum(preFP==1 & FP[i,]==1)
    a10 <- sum(preFP==1 & FP[i,]==0)
    a01 <- sum(preFP==0 & FP[i,]==1)
    
    R[i] <- a11/(a11+a10+a01)
    
  }
  
  # --------------------
  #   係数の値の重複
  # --------------------
  R1 <- unique(R)
  
  Rindex <- NULL
  for(i in 1:length(R1)){
    Rindex[i] <- list(which(R==R1[i])  )
  }
  
  G <- NULL
  for(i in 1:length(R1)){
    g <- 1
    if (length(Rindex[[i]])>1){
      fp <- FP[Rindex[[i]],] # 
      U1 <- unique(fp) # グループIDを付ける 
      for(k in 1:nrow(U1)){
        for(l in 1:nrow(fp)){
          if (sum(fp[l,]==U1[k,])==881){
            g[l] <- k
          }
        }
      }
    }
    G[i] <- list(g)
  }
  
  # --------------------
  #   ランキング結果
  # --------------------
  r.index <- order(R1,decreasing=TRUE)
  
  R2 <- R1[r.index]
  G2 <- G[r.index]
  
  L_CID <- NULL
  for(i in 1:length(Rindex)){
    L_CID[[i]] <- list(CID[Rindex[[i]]]) # CID
  }
  
  L_CID2 <- L_CID[r.index]
  
  # ---------------
  #   最終結果
  # ---------------
  ALL <- list(R2,G2,L_CID2) # Tanimoto係数,グループindex, CID
  
  # // 返り値
  return(ALL)
  
}

# -------------------------
# 	 ダウンロード無しver
# -------------------------
func_pub2b <- function(preID){
  
  library(rcdk)

  # ダウンロード
  url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",preID,"/SDF",sep="")
  download.file(url, "C:/R/pub.sdf", method = "auto", quiet=FALSE)
  mols <- load.molecules('C:/R/pub.sdf')
  fp <- get.fingerprint(mols[[1]], type="pubchem")
  
  preFP <- as.numeric(matrix(0,881)) # 初期値
  preFP[fp@bits] <- 1
  
  # -------------------------------
  #   sdfファイルのダウンロード
  # -------------------------------
  # PubChemで分子式検索
  x <- read.csv("C:/R/pub1.txt", header=FALSE)
  CID <- x[[1]]
  
  FP <- NULL
  for(i in 1:nrow(x)){
    mols <- load.molecules(paste0("C:/R/sdf/",CID[i],".sdf"))
    fp <- get.fingerprint(mols[[1]], type="pubchem")
    
    z <- as.numeric(matrix(0,881)) # 初期値
    z[fp@bits] <- 1
    
    FP <- rbind(FP,z)
  }
  
  # -----------------
  #   Tanimoto係数
  # -----------------
  R <- NULL
  for(i in 1:nrow(FP)){
    
    a11 <- sum(preFP==1 & FP[i,]==1)
    a10 <- sum(preFP==1 & FP[i,]==0)
    a01 <- sum(preFP==0 & FP[i,]==1)
    
    R[i] <- a11/(a11+a10+a01)
    
  }
  
  # --------------------
  #   係数の値の重複
  # --------------------
  R1 <- unique(R)
  
  Rindex <- NULL
  for(i in 1:length(R1)){
    Rindex[i] <- list(which(R==R1[i])  )
  }
  
  G <- NULL
  for(i in 1:length(R1)){
    g <- 1
    if (length(Rindex[[i]])>1){
      fp <- FP[Rindex[[i]],] # 
      U1 <- unique(fp) # グループIDを付ける 
      for(k in 1:nrow(U1)){
        for(l in 1:nrow(fp)){
          if (sum(fp[l,]==U1[k,])==881){
            g[l] <- k
          }
        }
      }
    }
    G[i] <- list(g)
  }
  
  # --------------------
  #   ランキング結果
  # --------------------
  r.index <- order(R1,decreasing=TRUE)
  
  R2 <- R1[r.index]
  G2 <- G[r.index]
  
  L_CID <- NULL
  for(i in 1:length(Rindex)){
    L_CID[[i]] <- list(CID[Rindex[[i]]]) # CID
  }
  
  L_CID2 <- L_CID[r.index]
  
  # ---------------
  #   最終結果
  # ---------------
  ALL <- list(R2,G2,L_CID2) # Tanimoto係数,グループindex, CID
  
  # // 返り値
  return(ALL)
  
}

# -----------------
# 	PubChem検索
# -----------------
func_pub2c <- function(CID){

  # ---------------------
  #   ファイル一覧削除
  # ---------------------
  foldername <- "C:/R/sdfrank/"
  flist <- paste0(foldername,list.files(foldername))
  
  if (length(list.files(foldername))>0){
    tt <- file.remove(flist)  
  }
  
  # ----------------------
  #   SDFファイル保存
  # ----------------------
  for(i in 1:length(CID)){
    url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",CID[i],"/SDF",sep="")
    download.file(url, paste0("C:/R/sdfrank/", CID[i], ".sdf"), method = "auto", quiet=FALSE)
  }

}