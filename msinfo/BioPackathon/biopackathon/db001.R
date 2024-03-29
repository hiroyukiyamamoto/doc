rm(list=ls(all=TRUE))

library(metabolomicsWorkbenchR)

# Introduction_to_metabolomicsWorkbenchR
# https://bioconductor.org/packages/devel/bioc/vignettes/metabolomicsWorkbenchR/inst/doc/Introduction_to_metabolomicsWorkbenchR.html

# 入力：試験タイトルをCancerで検索
# 出力：summary

df = do_query(
  context = 'study',
  input_item = 'study_title',
  input_value = 'Cancer',
  output_item = 'summary'
)

df = do_query(
  context = 'study',
  input_item = 'study_title',
  input_value = 'Schizophrenia',
  output_item = 'summary'
)



# Metabolomics Workbench REST URL-based API Specification
#  context : “study”, “compound”, “refmet”, “gene”, “protein”, “moverz” and “exactmass”

# Study IDがSTから始まるものを取得(全試験)
df = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST',
  output_item = 'summary'
)

# 全試験のStudy ID
allstudy_ids <- df$study_id
length(allstudy_ids) # 1256試験


# Project ID:PRから始まる
# Study ID:STから始まる

# 1つのProjectに対して、複数のStudy IDがついている


##  各種集計
#   生物種
pie(sort(table(df$subject_species)),radius=2)

