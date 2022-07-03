rm(list=ls(all=TRUE))

library(metabolomicsWorkbenchR)

# Introduction_to_metabolomicsWorkbenchR
# https://bioconductor.org/packages/devel/bioc/vignettes/metabolomicsWorkbenchR/inst/doc/Introduction_to_metabolomicsWorkbenchR.html

# Metabolomics Workbench REST URL-based API Specification
#  context : gstudyh, gcompoundh, grefmeth, ggeneh, gproteinh, gmoverzh and gexactmassh

# Study ID‚ªST‚©‚çn‚Ü‚é‚à‚Ì‚ğæ“¾(‘SŒ±)
df = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST',
  output_item = 'summary'
)

# ‘SŒ±‚ÌStudy ID
allstudy_ids <- df$study_id
length(allstudy_ids) # 1256Œ±

# Project ID:PR‚©‚çn‚Ü‚é
# Study ID:ST‚©‚çn‚Ü‚é

##  ŠeíWŒv
#   ¶•¨í
# pie(sort(table(df$subject_species)),radius=2)

# --------------------
#   ŒÂ•Êƒf[ƒ^‚Ìæ“¾
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
#   ‘Sƒf[ƒ^‚Ìæ“¾
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

# “®‚©‚È‚¢
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


