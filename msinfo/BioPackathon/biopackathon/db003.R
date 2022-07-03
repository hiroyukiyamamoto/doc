rm(list=ls(all=TRUE))

library(metabolomicsWorkbenchR)

# Introduction_to_metabolomicsWorkbenchR
# https://bioconductor.org/packages/devel/bioc/vignettes/metabolomicsWorkbenchR/inst/doc/Introduction_to_metabolomicsWorkbenchR.html

# --------------------------------------------
# "data" : “study_id” or “analysis_id”
# “datatable” : “analysis_id”
#  “mwtab” : “study_id” or “analysis_id”

# dataについては、そもそものAPIの仕様がstudy_idに複数のanalysis_idが付いているとエラーになる？
# mwtabでstudy_idで取得できない理由は謎
# analysis idでdatatableを取得してみる
# ---------------------------------------------
# ドキュメントをもう一度確認(APIとパッケージ両方)
# bioCのMLにメールする

# 全てのstudy IDを取得
# 各IDに対してanalysis IDを取得
# analysis IDからdatatableを取得する

# analysis idでは、検索ができない
# input_value = 'AN'ではエラーになる

# study idをすべて取得してから、analysis idを取得し、
# analysis idに対してdatatableを取得するのが良さそう

# step1 : Study IDがSTから始まるものを取得(全試験)
df = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST',
  output_item = 'summary'
)

# 全study IDを取得
allstudy_ids <- df$study_id

# -----------------------
#   DATA from study ID
# -----------------------
DF_data <- NULL
for(i in 1:length(allstudy_ids)){
  
  print(i)
  
  df = do_query(
    context = 'study',
    input_item = 'study_id',
    input_value = allstudy_ids[i],
    output_item = 'data'
  )
  
  DF_data[i] <- list(df)
}

# --------------------------
#   DATA from analysis ID
# --------------------------
# analysis IDを取得
all_ids <- NULL
for(i in 1:length(allstudy_ids)){
  
  df = do_query(
    context = 'study',
    input_item = 'study_id',
    input_value = allstudy_ids[i],
    output_item = 'analysis'
  )  
  
  all_ids <- rbind(all_ids,cbind(df$study_id,df$analysis_id))
  
}

allflag <- NULL
datatable_all <- NULL;mwtab_all <- NULL
for(i in 1:nrow(all_ids)){
  
  print(i)
  errflag <- c(0,0)
  
  df_data <- NULL
  flag <- tryCatch({
    flag <- 0
    df_data = do_query(
      context = 'study',
      input_item = 'analysis_id',
      input_value = all_ids[i,2],
      output_item = 'datatable')
    flag <- 0
    },
    error = function(e) {
      message("ERROR in datatable!")
      return(1)
    }
  )
  
  if (flag==1){
    errflag[1] <- 1
  }
  
  #datatable_all[i] <- list(df_data)
  
  # mwtabの場合
  mwtab_data <- NULL
  flag <- tryCatch({
    df_mwtab = do_query(
      context = 'study',
      input_item = 'analysis_id',
      input_value = all_ids[i,2],
      output_item = 'mwtab')
    flag <- 0
    },
    error=function(e){
      message("ERROR in mwtab!")
      return(1)
    }
  )
  
  if (flag==1){
    errflag[2] <- 1
  }

  
  print(errflag)
  
  #mwtab_all[i] <- list(df_mwtab)
  allflag <- rbind(allflag,c(all_ids[i,],errflag))
}

# ---------------------------------------------------------------------

table_index <- NULL; mwtab_index <- NULL
for(i in 1:length(datatable_all)){
  table_index[i] <- is.null(datatable_all[[i]])
  mwtab_index[i] <- is.null(mwtab_all[[i]])
}


# 各analysis IDに対してデータを取得する




df_data = do_query(
  context = 'study',
  input_item = 'analysis_id',
  input_value = 'AN000076',
  output_item = 'datatable'
)



# 全試験のStudy ID
allstudy_ids <- df$study_id
length(allstudy_ids) # 1256試験

# Project ID:PRから始まる
# Study ID:STから始まる

##  各種集計
#   生物種
# pie(sort(table(df$subject_species)),radius=2)

# --------------------
#   個別データの取得
# --------------------
df_data = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST000046',
  output_item = 'data'
)

df_mwtab = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST000001',
  output_item = 'mwtab'
)



# --------------------
#   全データの取得
# --------------------
DF_data <- NULL
for(i in 1:length(allstudy_ids)){
  
  print(i)
  
  df = do_query(
    context = 'study',
    input_item = 'study_id',
    input_value = allstudy_ids[i],
    output_item = 'data'
  )
  
  DF_data[i] <- list(df)
}

# 動かない
#df = do_query(
#  context = 'study',
#  input_item = 'study_id',
#  input_value = 'ST',
#  output_item = 'datatable'
#)

#DF_mwtab <- NULL
#for(i in 1:length(allstudy_ids)){
#for(i in 1:100){  
#  df = do_query(
#    context = 'study',
#    input_item = 'study_id',
#    input_value = allstudy_ids[i],
#    output_item = 'mwtab'
#  )
#  
#  DF_mwtab[i] <- list(df)
#}

# save(DF_data,file="C:/Users/ito/Documents/R/biopackathon/DF.RData")

index <- NULL
for(i in 1:length(DF_data)){
  index[i] <- is.null(DF_data[[i]])
}


AN000076




df_data = do_query(
  context = 'study',
  input_item = 'analysis_id',
  input_value = all_ids[i,2],
  output_item = 'mwtab'
)

df_data = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = all_ids[i,1],
  output_item = 'data'
)

