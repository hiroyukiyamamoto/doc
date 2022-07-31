### リピドミクス

- 解析例　LC-MS/MS (DDA)
  - https://hiroyukiyamamoto.github.io/doc/msinfo/msinfo001.html
  - 必要なライブラリ
    - mzR : .mzmlファイルの読み込み
    - MsBackendMsp : .mspファイルの読み込み
    - MSnbase : MS/MSスペクトルの処理
  - 解析の流れ
    - [事前準備] MS1のデータからピークを拾い、ターゲットとなるピークのm/zとRTを取得する
    - ターゲットピークのMS/MSスペクトルを抽出する (ターゲットピークのm/zと、precursorのm/zの一致を調べる)
    - .mspから該当するスペクトルライブラリを抽出する (precursor m/zの情報から、)
