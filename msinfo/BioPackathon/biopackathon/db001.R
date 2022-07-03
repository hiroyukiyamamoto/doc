rm(list=ls(all=TRUE))

library(metabolomicsWorkbenchR)

# Introduction_to_metabolomicsWorkbenchR
# https://bioconductor.org/packages/devel/bioc/vignettes/metabolomicsWorkbenchR/inst/doc/Introduction_to_metabolomicsWorkbenchR.html

# “ü—ÍFŒ±ƒ^ƒCƒgƒ‹‚ğCancer‚ÅŒŸõ
# o—ÍFsummary

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

# 1‚Â‚ÌProject‚É‘Î‚µ‚ÄA•¡”‚ÌStudy ID‚ª‚Â‚¢‚Ä‚¢‚é


##  ŠeíWŒv
#   ¶•¨í
pie(sort(table(df$subject_species)),radius=2)

