### プロテオミクス

- 解析例　
  - 検討中
  - 必要なライブラリ
    - mzR : .mzmlファイルの読み込み
    - MSnbase : MS/MSスペクトルの処理
  - 解析手順
    - 事前準備
      - 検討中
    - 検討中
- (メモ) 解析用データ　
  - https://repository.jpostdb.org/entry/JPST001257
  - rawファイル、msfファイルがある
    - F150909_FFPE_Nikkyo_180min_MS1_HCDIT_MS2_HCDOT_5uL_Sample01.raw
    - F150909_FFPE_Nikkyo_180min_MS1_HCDIT_MS2_HCDOT_5uL_Sample01.msf
    - DDA,トリプシン処理
      - msf : Proteome discovererの解析結果のファイル
      - msfをmgfに変換してmascotに投げる
  - https://bioconductor.org/packages/release/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.html
  - https://bioconductor.riken.jp/packages/3.0/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.pdf
  - 以前はrTandemが使われていたが、現在はBioconductorから削除されている
  - 現在は、ソフトウェアもしくはデータベース(Mascot, Comet, MaxQuant, MSGF+, OMSSAなど)で検索し、その結果をRで読み込んで利用するのが標準的
  - While searches are generally performed using third-party software independently of R or can be started from R using a system call, the rTANDEM package allows one to execute such searches using the X!Tandem engine. 
    - https://bioconductor.riken.jp/packages/3.8/workflows/vignettes/proteomics/inst/doc/proteomics.html#msms-database-search
  - 検討事項
    - ソフトウェアで一通り解析し、その結果をRで読み込む
    - MS/MSスペクトルをmgfファイルに出力して、Mascotなどのデータベースで検索する。検索結果を再度Rで読み込む
- 解析用メモ
  - 特定のペプチドに着目してDDAのスペクトルを取得し、mgfファイルにしてMascotに投げる
  - msfファイルを確認し、注目するペプチドを決めてMascotに投げる   

  
