# -----------------
# 	PubChem����
# -----------------
func_pub2 <- function(mf,preID){

  library(rcdk)
  
  # ---------------------
  #   �t�@�C���ꗗ�폜
  # ---------------------
  foldername <- "C:/R/sdf/"
  flist <- paste0(foldername,list.files(foldername))
  
  if (length(list.files(foldername))>0){
    tt <- file.remove(flist)  
  }
  
  # -----------------------
  #   ���q������\������
  # -----------------------
  # mf <- "C7H12N2O5" # ���m�s�[�N�̕��q��
  
  # PubChem�ŕ��q������
  url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/",mf,"/cids/TXT",sep="")
  
  # CID�̃��X�g���_�E�����[�h
  download.file(url, "C:/R/pub1.txt", method = "auto", quiet=FALSE)
  
  # -------------------------------
  #   sdf�t�@�C���̃_�E�����[�h
  # -------------------------------
  # PubChem�ŕ��q������
  x <- read.csv("C:/R/pub1.txt", header=FALSE)
  CID <- x[[1]]
  
  FP <- NULL
  for(i in 1:nrow(x)){
    url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",CID[i],"/SDF",sep="")
    download.file(url, paste0("C:/R/sdf/", CID[i], ".sdf"), method = "auto", quiet=FALSE)
    mols <- load.molecules(paste0("C:/R/sdf/",CID[i],".sdf"))
    fp <- get.fingerprint(mols[[1]], type="pubchem")
    
    z <- as.numeric(matrix(0,881)) # �����l
    z[fp@bits] <- 1
    
    FP <- rbind(FP,z)
  }
  
  # -------------------------------
  #   �O��̂�pubchem fingerprint
  # -------------------------------
  # pubchem CID�i���������X�g�j
  # preID <- 119 # �v�ݒ�
  
  # �_�E�����[�h
  url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",preID,"/SDF",sep="")
  download.file(url, "C:/R/pub.sdf", method = "auto", quiet=FALSE)
  mols <- load.molecules('C:/R/pub.sdf')
  fp <- get.fingerprint(mols[[1]], type="pubchem")
  
  preFP <- as.numeric(matrix(0,881)) # �����l
  preFP[fp@bits] <- 1
  
  # -----------------
  #   Tanimoto�W��
  # -----------------
  R <- NULL
  for(i in 1:nrow(FP)){
    
    a11 <- sum(preFP==1 & FP[i,]==1)
    a10 <- sum(preFP==1 & FP[i,]==0)
    a01 <- sum(preFP==0 & FP[i,]==1)
    
    R[i] <- a11/(a11+a10+a01)
    
  }
  
  # --------------------
  #   �W���̒l�̏d��
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
      U1 <- unique(fp) # �O���[�vID��t���� 
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
  #   �����L���O����
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
  #   �ŏI����
  # ---------------
  ALL <- list(R2,G2,L_CID2) # Tanimoto�W��,�O���[�vindex, CID
  
  # // �Ԃ�l
  return(ALL)
  
}

# -------------------------
# 	 �_�E�����[�h����ver
# -------------------------
func_pub2b <- function(preID){
  
  library(rcdk)

  # �_�E�����[�h
  url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",preID,"/SDF",sep="")
  download.file(url, "C:/R/pub.sdf", method = "auto", quiet=FALSE)
  mols <- load.molecules('C:/R/pub.sdf')
  fp <- get.fingerprint(mols[[1]], type="pubchem")
  
  preFP <- as.numeric(matrix(0,881)) # �����l
  preFP[fp@bits] <- 1
  
  # -------------------------------
  #   sdf�t�@�C���̃_�E�����[�h
  # -------------------------------
  # PubChem�ŕ��q������
  x <- read.csv("C:/R/pub1.txt", header=FALSE)
  CID <- x[[1]]
  
  FP <- NULL
  for(i in 1:nrow(x)){
    mols <- load.molecules(paste0("C:/R/sdf/",CID[i],".sdf"))
    fp <- get.fingerprint(mols[[1]], type="pubchem")
    
    z <- as.numeric(matrix(0,881)) # �����l
    z[fp@bits] <- 1
    
    FP <- rbind(FP,z)
  }
  
  # -----------------
  #   Tanimoto�W��
  # -----------------
  R <- NULL
  for(i in 1:nrow(FP)){
    
    a11 <- sum(preFP==1 & FP[i,]==1)
    a10 <- sum(preFP==1 & FP[i,]==0)
    a01 <- sum(preFP==0 & FP[i,]==1)
    
    R[i] <- a11/(a11+a10+a01)
    
  }
  
  # --------------------
  #   �W���̒l�̏d��
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
      U1 <- unique(fp) # �O���[�vID��t���� 
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
  #   �����L���O����
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
  #   �ŏI����
  # ---------------
  ALL <- list(R2,G2,L_CID2) # Tanimoto�W��,�O���[�vindex, CID
  
  # // �Ԃ�l
  return(ALL)
  
}

# -----------------
# 	PubChem����
# -----------------
func_pub2c <- function(CID){

  # ---------------------
  #   �t�@�C���ꗗ�폜
  # ---------------------
  foldername <- "C:/R/sdfrank/"
  flist <- paste0(foldername,list.files(foldername))
  
  if (length(list.files(foldername))>0){
    tt <- file.remove(flist)  
  }
  
  # ----------------------
  #   SDF�t�@�C���ۑ�
  # ----------------------
  for(i in 1:length(CID)){
    url <- paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",CID[i],"/SDF",sep="")
    download.file(url, paste0("C:/R/sdfrank/", CID[i], ".sdf"), method = "auto", quiet=FALSE)
  }

}