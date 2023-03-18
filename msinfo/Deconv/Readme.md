### MS/MSデコンボリューション

- osdパッケージを用いたデコンボリューション
  - https://hiroyukiyamamoto.github.io/doc/msinfo/Deconv/deconv001.html
  - https://cran.r-project.org/web/packages/osd/index.html
- ALSパッケージ
  - 確認中
  - https://cran.r-project.org/web/packages/ALS/index.html 
- alsaceパッケージ
  - 要確認
  - http://bioconductor.jp/packages/3.11/bioc/html/alsace.html
- 解析手法(実装予定)
  - デコンボリューション
    - [済]ALS
    - [済]RALS(Regularized ALS)、ACLS
    - [済]non-negative ALS
    - ICA
    - nonnegative ICA
    - non-negative matrix factorization
  - 制約条件
    - [済]normalize
    - unimodal
  - 成分数の推定
    - singular value ratio (SVR) 
    - IND, RESO
  - データ
    - HILIC-Pos-SWATH-25Da-20140701_08_GB004467_Swath25Da.mzML
- 参考資料
  - 論文
    - Hiroyuki Yamamoto, Keishi Hada, Hideki Yamaji, Tomohisa Katsuda, Hiromu Ohno, Hideki Fukuda,
    - "Application of regularized alternating least squares and independent component analysis to HPLC-DAD data of Haematococcus pluvialis metabolites", Biochem. Eng. Journal, 32 (2006) 149-156.
    - https://www.sciencedirect.com/science/article/abs/pii/S1369703X06002610?via%3Dihub
  - 学会発表
    - ALS と ICA のスペクトル分離法への応用 Haematococcus pluvialis代謝物質のHPLC-DADデータの解析
    - 第49回自動制御連合講演会 2006 年 11 月 25 日，26 日 神戸大学
    - https://www.jstage.jst.go.jp/article/jacc/49/0/49_0_389/_article/-char/ja/
  - 参考文献
    - Chen, Zeng‐Ping et al. “Determination of the number of components in mixtures using a new approach incorporating chemical information.” Journal of Chemometrics 13 (1999): n. pag.
    - https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/pdf/10.1002/%28SICI%291099-128X%28199901/02%2913%3A1%3C15%3A%3AAID-CEM527%3E3.0.CO%3B2-I






