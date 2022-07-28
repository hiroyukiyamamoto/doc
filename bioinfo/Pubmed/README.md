### Pubmedデータ取得(R)

- Rパッケージの調査
    - easyPubMed 
    - RISmed
- 雑感(easyPubMed)
    - easyPubmedのget_pubmed_ids関数は、PubmedIDを20個までしか取得できない
    - fetch_all_pubmed_ids関数で取得しようとしても、全データを取得できない(件数が多いせいか？)
    - Note that only 0 PubMed IDs were retrievedメッセージ
- 雑感(RISmed)
    - 高速で簡単に取得できる
    - プログラム(例)
        - library(RISmed)
        - keyword <- '"neutral fat" OR "triglyceride"'
        - all <- EUtilsSummary(keyword,type="esearch",db="pubmed")
        - counts <- QueryCount(all)
        - all <- EUtilsSummary(keyword,type="esearch",db="pubmed", retmax=counts)
        - pmid <- QueryId(all)
        - which(pmid=="25565485")


