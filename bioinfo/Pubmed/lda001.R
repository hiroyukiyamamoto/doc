

# https://rpubs.com/gingi99/16714

# ベクトル型で以下のように格納できるデータならOK。csvファイルなら結合していけばOK。
sentence <- c("I am the very model of a modern major general", "I have a major headache")
# lexicalize関数で、LDAで分析するためのデータを生成。LIST型として展開される。
# LISTの要素は2つの成分。LIST型の$documentsと文字列ベクトルの$vovabに。
test <- lexicalize(sentence, lower = TRUE)
# testの中身を見ると、$documentはセンテンスの出現位置と回数、$vocabがユニークな単語
test


data(cora.documents)
head(cora.documents, n = 2)

data(cora.vocab)
head(cora.vocab)

data(cora.titles)
head(cora.titles)

# 分析データの作成（トリッキーな参照をしているので注意）
# 1列目がcora.documentsの第一成分で使われる単語のリスト、2列目がその出現回数
data_cora <- as.data.frame(cbind(cora.vocab[cora.documents[[1]][1, ] + 1], cora.documents[[1]][2, 
]))
# coreの1番目の記事はこれらの単語とその出現回数で構成されていることが分かる。
head(data_cora)

# 推定するトピック数の設定
k <- 10

# ldaパッケージはギブスサンプラーでやるようです。
# ギブスサンプラーの中でも3つくらいmethodがあるようです。
result <- lda.collapsed.gibbs.sampler(cora.documents, 
                                      k,
                                      cora.vocab,
                                      25,  # 繰り返し数
                                      0.1, # ディリクレ過程のハイパーパラメータα
                                      0.1, # ディリクレ過程のハイパーパラメータη
                                      compute.log.likelihood=TRUE)

# サマリを見ると、10成分のリストで構成されている。
# assignments：文書Dと同じ長さのリスト。値は単語が割り当てられたトピックNoを示す。
# topic：k × vの行列。値はそのトピックに出現する単語数を表す。
# topic_sums：それぞれのトピックに割り当てられた単語の合計数
# document_sums：k × Dの行列。割り振られたトピックにおける一文章内の単語数を示す。
summary(result)




