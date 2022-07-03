rm(list=ls(all=TRUE))

library(metabolomicsWorkbenchR)

# Introduction_to_metabolomicsWorkbenchR
# https://bioconductor.org/packages/devel/bioc/vignettes/metabolomicsWorkbenchR/inst/doc/Introduction_to_metabolomicsWorkbenchR.html

# ���́F�����^�C�g����Cancer�Ō���
# �o�́Fsummary

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
#  context : �gstudy�h, �gcompound�h, �grefmet�h, �ggene�h, �gprotein�h, �gmoverz�h and �gexactmass�h

# Study ID��ST����n�܂���̂��擾(�S����)
df = do_query(
  context = 'study',
  input_item = 'study_id',
  input_value = 'ST',
  output_item = 'summary'
)

# �S������Study ID
allstudy_ids <- df$study_id
length(allstudy_ids) # 1256����


# Project ID:PR����n�܂�
# Study ID:ST����n�܂�

# 1��Project�ɑ΂��āA������Study ID�����Ă���


##  �e��W�v
#   ������
pie(sort(table(df$subject_species)),radius=2)
