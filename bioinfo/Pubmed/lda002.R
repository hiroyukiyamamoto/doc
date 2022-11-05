

library(topicmodels)

data("AssociatedPress", package = "topicmodels")
lda <- LDA(AssociatedPress[1:20,], control = list(alpha = 0.1), k = 2)
lda_inf <- posterior(lda, AssociatedPress[21:30,])


# latent semantic analysis
# https://cran.r-project.org/web/packages/lsa/lsa.pdf
# https://cran.r-project.org/web/views/NaturalLanguageProcessing.html

# --------
#   LSA
# --------
library(lsa)
td = tempfile()
dir.create(td)
write( c("dog", "cat", "mouse"), file=paste(td, "D1", sep="/") )
write( c("ham", "mouse", "sushi"), file=paste(td, "D2", sep="/") )
write( c("dog", "pet", "pet"), file=paste(td, "D3", sep="/") )

data(stopwords_en)
myMatrix = textmatrix(td, stopwords=stopwords_en) # 単純なbag-of-the wordsのmatrixに変換？
myMatrix = lw_logtf(myMatrix) * gw_idf(myMatrix) # 重み付けがされている？
myLSAspace = lsa(myMatrix, dims=dimcalc_share()) # LSA
as.textmatrix(myLSAspace)
# clean up
unlink(td, recursive=TRUE)
