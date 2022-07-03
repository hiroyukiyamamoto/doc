rm(list=ls(all=TRUE))

library(metabolomicsWorkbenchR)

# Introduction_to_metabolomicsWorkbenchR
# https://bioconductor.org/packages/devel/bioc/vignettes/metabolomicsWorkbenchR/inst/doc/Introduction_to_metabolomicsWorkbenchR.html

# --------------------------------------------
# "data" : �gstudy_id�h or �ganalysis_id�h
# �gdatatable�h : �ganalysis_id�h
#  �gmwtab�h : �gstudy_id�h or �ganalysis_id�h

# data�ɂ��ẮA����������API�̎d�l��study_id�ɕ�����analysis_id���t���Ă���ƃG���[�ɂȂ�H
# mwtab��study_id�Ŏ擾�ł��Ȃ����R�͓�
# analysis id��datatable���擾���Ă݂�
# ---------------------------------------------
# �h�L�������g��������x�m�F(API�ƃp�b�P�[�W����)
# bioC��ML�Ƀ��[������

# �S�Ă�study ID���擾
# �eID�ɑ΂���analysis ID���擾
# analysis ID����datatable���擾����

# analysis id�ł́A�������ł��Ȃ�
# input_value = 'AN'�ł̓G���[�ɂȂ�

# study id�����ׂĎ擾���Ă���Aanalysis id���擾���A
# analysis id�ɑ΂���datatable���擾����̂��ǂ�����

# step1 : Study ID��ST����n�܂���̂��擾(�S����)
df = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST',
  output_item = 'summary'
)

# �Sstudy ID���擾
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
# analysis ID���擾
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
  
  # mwtab�̏ꍇ
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


# �eanalysis ID�ɑ΂��ăf�[�^���擾����




df_data = do_query(
  context = 'study',
  input_item = 'analysis_id',
  input_value = 'AN000076',
  output_item = 'datatable'
)



# �S������Study ID
allstudy_ids <- df$study_id
length(allstudy_ids) # 1256����

# Project ID:PR����n�܂�
# Study ID:ST����n�܂�

##  �e��W�v
#   ������
# pie(sort(table(df$subject_species)),radius=2)

# --------------------
#   �ʃf�[�^�̎擾
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
#   �S�f�[�^�̎擾
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

# �����Ȃ�
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
