# doc2vec

library(doc2vec)
library(tokenizers.bpe)
library(udpipe)
data(belgium_parliament, package = "tokenizers.bpe")
x <- subset(belgium_parliament, language %in% "dutch")
x <- data.frame(doc_id = sprintf("doc_%s", 1:nrow(x)), 
                text   = x$text, 
                stringsAsFactors = FALSE)
x$text   <- tolower(x$text)
x$text   <- gsub("[^[:alpha:]]", " ", x$text)
x$text   <- gsub("[[:space:]]+", " ", x$text)
x$text   <- trimws(x$text)
x$nwords <- txt_count(x$text, pattern = " ")
x        <- subset(x, nwords < 1000 & nchar(text) > 0)

# --------------------
#   Build the model
# --------------------
## Low-dimensional model using DM, low number of iterations, for speed and display purposes
model <- paragraph2vec(x = x, type = "PV-DM", dim = 5, iter = 3,  
                       min_count = 5, lr = 0.05, threads = 1)
str(model)


## More realistic model
model <- paragraph2vec(x = x, type = "PV-DBOW", dim = 100, iter = 20, 
                       min_count = 5, lr = 0.05, threads = 4)

# --------------------------------------------------------------------------------
#   Get similar documents or words when providing sentences, documents or words
# --------------------------------------------------------------------------------
sentences <- strsplit(setNames(x$text, x$doc_id), split = " ")
nn <- predict(model, newdata = sentences, type = "nearest", which = "sent2doc", top_n = 5)

# ------------
#   top2vec
# ------------
library(doc2vec)
library(word2vec)
library(uwot)
library(dbscan)
data(be_parliament_2020, package = "doc2vec")
x      <- data.frame(doc_id = be_parliament_2020$doc_id,
                     text   = be_parliament_2020$text_nl,
                     stringsAsFactors = FALSE)
x$text <- txt_clean_word2vec(x$text)
x      <- subset(x, txt_count_words(text) < 1000)

d2v    <- paragraph2vec(x, type = "PV-DBOW", dim = 50, 
                        lr = 0.05, iter = 10,
                        window = 15, hs = TRUE, negative = 0,
                        sample = 0.00001, min_count = 5, 
                        threads = 1)

#str(model)
#embedding <- as.matrix(model, which = "words")
#embedding <- as.matrix(model, which = "docs")
#head(embedding)

#vocab <- summary(model, type = "vocabulary",  which = "docs")
#vocab <- summary(model, type = "vocabulary",  which = "words")


# d2vに対してpredictすればよい
# umapを見るときはtop2vec

model  <- top2vec(d2v, 
                  control.dbscan = list(minPts = 50), 
                  control.umap = list(n_neighbors = 15L, n_components = 3), umap = tumap, 
                  trace = TRUE)

plot(model$umap, pch=16, cex=0.1)
plot(model$umap, pch=16, cex=0.1, col=model$dbscan$cluster+1)


info   <- summary(model, top_n = 7)
info$topwords
