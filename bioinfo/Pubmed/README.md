### Pubmedデータ取得(R)

- Rパッケージの調査
    - easyPubMed 
    - RISmed
- 雑感(easyPubMed)　いまいち
    - easyPubmedのget_pubmed_ids関数は、PubmedIDを20個までしか取得できない
    - fetch_all_pubmed_ids関数で取得しようとしても、全データを取得できない(件数が多いせいか？)
    - Note that only 0 PubMed IDs were retrievedメッセージ
- 雑感(RISmed)　とても良い
    - 高速で簡単に取得できる
    - プログラム(例)
        - library(RISmed)
        - keyword <- '"neutral fat" OR "triglyceride"'
        - all <- EUtilsSummary(keyword,type="esearch",db="pubmed")
        - counts <- QueryCount(all) # 件数を先に取得
        - all <- EUtilsSummary(keyword,type="esearch",db="pubmed", retmax=counts) # 件数を指定して取得
        - pmid <- QueryId(all)
        - which(pmid=="25565485")
   - EUtilsSummaryのretmaxは大きめに設定するとエラーになるので、具体的な件数を指定する
- メモ
    - 要検討 fobitools
    - fooDB、食レポ
    - PubMedWordcloud
    - 特定の現象で検索、特定の代謝物で検索、類似度データを生成、解析
 


