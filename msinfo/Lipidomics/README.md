# 前提と準備

## R パッケージ

以下の Bioconductor
パッケージが必要です。インストール済みの場合は省略できます。

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "xcms",
  "CAMERA",
  "MSnbase",
  "Spectra",
  "CompoundDb",
  "S4Vectors"
))
```

## 外部データ

各 script の先頭の設定ブロックでパスを指定します。

-   **mzML ファイル**: LC-MS/MS の生データを mzML 形式に変換したもの
-   **LIPID MAPS LMSD**: [LIPID
    MAPS](https://www.lipidmaps.org/databases/lmsd) からダウンロードした
    `structures.sdf`
-   **LipidBlast ライブラリ**: MSP 形式の MS2 スペクトルライブラリ

# ワークフローの全体像

この解析は、大きく次の 9 段階に分かれます。

1.  `mzML` からサンプル比較用のピーク候補表を作る
2.  同じ化合物に由来しそうなピークをまとめる
3.  脂質候補として優先して見るピークを絞る
4.  `LIPID MAPS` を使って MS1 の質量から候補名を付ける
5.  `LipidBlast` から必要なライブラリ候補だけを取り出す
6.  実測 DDA MS2 をピーク候補ごとにまとめる
7.  `compareSpectra()` で実測 MS2 とライブラリを比べる
8.  MS1 と MS2 の結果を 1 つの annotation 表にまとめる
9.  PCA と loadings でサンプル差と寄与脂質を読む

# スクリプト一覧とコード

## 1. XCMS でピーク候補の表を作る

mzML の生データからピーク検出・アライメントを行い、
サンプル間で比較できる feature 表と intensity 行列を作ります。

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# このスクリプトでは、複数の mzML ファイルを XCMS で一括処理し、
# feature 行列を作成した後に PCA を行う。
# さらに、PC1 と PC2 のスコアプロットを PNG 画像として保存する。

# xcms が利用できるか最初に確認する
if (!requireNamespace("xcms", quietly = TRUE)) {
  stop("Package 'xcms' is required.")
}

# 入力 mzML フォルダ
mzml_dir <- "./mzml"

# 結果を書き出すフォルダ
output_dir <- "./output/xcms_pca"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# mzML ファイル一覧を取得する
mzfiles <- list.files(mzml_dir, pattern = "\\.mzML$", full.names = TRUE)
if (length(mzfiles) == 0) {
  stop("No mzML files were found.")
}

# ファイル名から拡張子を除いたものをサンプル名として使う
sample_names <- tools::file_path_sans_ext(basename(mzfiles))
names(mzfiles) <- sample_names

# サンプル名から藻類タイプを取り出す
# 例:
#   Posi_Ida_Chlamydomonas_1      -> Chlamydomonas
#   Posi_Ida_Chrollera_Utex2341_1 -> Chrollera_Utex2341
extract_algae_type <- function(sample_name) {
  sub("_[0-9]+$", "", sub("^Posi_Ida_", "", sample_name))
}
```

### centWave でピーク検出を行う

``` r
# 1. centWave でピーク検出を行う
print(paste("Found", length(mzfiles), "mzML files"))
print("Running XCMS peak detection with centWave...")
xset <- xcms::xcmsSet(
  files = mzfiles,
  method = "centWave",
  ppm = 10,
  peakwidth = c(5, 20),
  snthresh = 3
)
```

### 保持時間補正の前に、一度ピークをグループ化する

``` r
# 2. 保持時間補正の前に、一度ピークをグループ化する
# ここでは各サンプル間で「同じ feature らしいピーク」をまとめる
print("Grouping peaks before retention time correction...")
xset <- xcms::group(
  xset,
  bw = 5,
  mzwid = 0.015,
  minfrac = 0.5
)
```

### 複数ファイルがある場合は保持時間補正を行う

``` r
# 3. 複数ファイルがある場合は保持時間補正を行う
# 補正後に再度 group() を実行して、feature 対応を整え直す
if (length(mzfiles) >= 2) {
  print("Performing retention time correction with obiwarp...")
  xset <- xcms::retcor(xset, method = "obiwarp")

  print("Grouping peaks after retention time correction...")
  xset <- xcms::group(
    xset,
    bw = 5,
    mzwid = 0.015,
    minfrac = 0.5
  )
}
```

### fillPeaks() で欠損ピークを補完する

``` r
# 4. fillPeaks() で欠損ピークを補完する
# これにより、サンプル間で同じ feature を比較しやすくなる
print("Filling missing peaks...")
xset <- xcms::fillPeaks(xset)
```

### XCMS の結果から intensity 行列を作る

``` r
# 5. XCMS の結果から intensity 行列を作る
# groupval() は「feature x sample」で返すため、
# PCA しやすいように転置して「sample x feature」に直す
print("Extracting integrated peak intensity matrix...")
feature_matrix <- xcms::groupval(
  xset,
  value = "into",
  method = "medret",
  intensity = "into"
)
feature_matrix <- t(as.matrix(feature_matrix))
feature_matrix[is.na(feature_matrix)] <- 0
rownames(feature_matrix) <- sample_names

# 各 feature の m/z や保持時間などの情報も別に保存する
feature_definitions <- xcms::groups(xset)
feature_ids <- paste0("feature_", seq_len(nrow(feature_definitions)))
colnames(feature_matrix) <- feature_ids
rownames(feature_definitions) <- feature_ids

print(
  paste(
    "Feature matrix dimensions:",
    nrow(feature_matrix), "samples x", ncol(feature_matrix), "features"
  )
)
```

### 全サンプルで強度が 0 の feature は削除する

``` r
# 6. 全サンプルで強度が 0 の feature は削除する
# こうした列は PCA に寄与しないため除いておく
nonzero_cols <- which(colSums(feature_matrix) > 0)
feature_matrix_filtered <- feature_matrix[, nonzero_cols, drop = FALSE]
feature_definitions_filtered <- feature_definitions[nonzero_cols, , drop = FALSE]
print(
  paste(
    "Retained", ncol(feature_matrix_filtered),
    "features after removing all-zero columns"
  )
)
```

### 強度差を少し穏やかにするために log1p 変換を行う

``` r
# 7. 強度差を少し穏やかにするために log1p 変換を行う
# log1p(x) = log(1 + x) なので、0 を含む行列でも扱いやすい
log_matrix <- log1p(feature_matrix_filtered)

pca_result <- NULL
pca_scores <- NULL
```

### PCA を計算する

``` r
# 8. PCA を計算する
# PCA はサンプル数・特徴量数ともに最低 2 以上必要
if (nrow(log_matrix) >= 2 && ncol(log_matrix) >= 2) {
  # 分散が 0 の列は PCA に不要なので除く
  nonzero_sd_cols <- which(apply(log_matrix, 2, sd) > 0)
  pca_input <- log_matrix[, nonzero_sd_cols, drop = FALSE]
  feature_definitions_pca <- feature_definitions_filtered[nonzero_sd_cols, , drop = FALSE]

  if (ncol(pca_input) >= 2) {
    print("Running PCA...")
    pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
    pca_scores <- as.data.frame(pca_result$x)
    pca_scores$sample <- rownames(pca_scores)
    pca_scores$algae_type <- extract_algae_type(pca_scores$sample)

    # PCA に使われた feature 定義も別保存しておく
    save(
      feature_definitions_pca,
      file = file.path(output_dir, "xcms_feature_definitions_for_pca.rds")
    )
  } else {
    warning("Too few variable features remained after filtering. PCA was skipped.")
  }
} else {
  warning("At least two mzML files are needed for PCA. Feature matrix was created, but PCA was skipped.")
}
```

### 中間結果と最終結果を保存する

``` r
# 9. 中間結果と最終結果を保存する
save(xset, file = file.path(output_dir, "xcms_object.rds"))
save(feature_matrix, file = file.path(output_dir, "xcms_feature_matrix_raw.rds"))
save(feature_matrix_filtered, file = file.path(output_dir, "xcms_feature_matrix_filtered.rds"))
save(feature_definitions_filtered, file = file.path(output_dir, "xcms_feature_definitions.rds"))

write.csv(
  feature_matrix,
  file = file.path(output_dir, "xcms_feature_matrix_raw.csv"),
  row.names = TRUE
)
write.csv(
  feature_matrix_filtered,
  file = file.path(output_dir, "xcms_feature_matrix_filtered.csv"),
  row.names = TRUE
)
write.csv(
  feature_definitions_filtered,
  file = file.path(output_dir, "xcms_feature_definitions.csv"),
  row.names = TRUE
)
```

### PCA が計算できた場合は、スコア表とスコアプロットも保存する

``` r
# 10. PCA が計算できた場合は、スコア表とスコアプロットも保存する
if (!is.null(pca_result)) {
  save(pca_result, file = file.path(output_dir, "xcms_pca_result.rds"))
  save(pca_scores, file = file.path(output_dir, "xcms_pca_scores.rds"))
  write.csv(
    pca_scores,
    file = file.path(output_dir, "xcms_pca_scores.csv"),
    row.names = FALSE
  )

  # 藻類タイプごとに色を割り当てる
  algae_types <- unique(pca_scores$algae_type)
  palette_colors <- grDevices::hcl.colors(length(algae_types), palette = "Dark 3")
  names(palette_colors) <- algae_types

  # PC1, PC2 の寄与率を軸ラベルに表示する
  pc1_var <- round(100 * summary(pca_result)$importance[2, 1], 1)
  pc2_var <- round(100 * summary(pca_result)$importance[2, 2], 1)

  # PC1-PC2 スコアプロットを PNG 保存する
  png(
    filename = file.path(output_dir, "xcms_pca_score_plot.png"),
    width = 1600,
    height = 1200,
    res = 150
  )
  plot(
    pca_scores$PC1,
    pca_scores$PC2,
    col = palette_colors[pca_scores$algae_type],
    pch = 19,
    cex = 1.3,
    xlab = paste0("PC1 (", pc1_var, "%)"),
    ylab = paste0("PC2 (", pc2_var, "%)"),
    main = "PCA Score Plot Colored by Algae Type"
  )
  text(
    pca_scores$PC1,
    pca_scores$PC2,
    labels = pca_scores$sample,
    pos = 3,
    cex = 0.6
  )
  legend(
    "topright",
    legend = algae_types,
    col = palette_colors[algae_types],
    pch = 19,
    bty = "n",
    cex = 0.9
  )
  dev.off()

  print("PCA completed and saved.")
}

# 最後に保存先を表示する
print(paste("Saved outputs to", output_dir))
```

## 2. CAMERA で同じ化合物に由来しそうなピークをまとめる

XCMS の結果をもとに、同じ化合物に由来する可能性が高い feature
をグループ化し、 同位体・アダクトの候補を付けます。

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# このスクリプトでは、mspp_xcms_pca.R で保存した xcms オブジェクトを読み込み、
# CAMERA を使って feature にアイソトープ・アダクト候補の注釈を付ける。
# まずは「化学式や化合物名の同定」ではなく、
# 「同じピークグループ内でどの feature が関連していそうか」を整理する段階を目的とする。

if (!requireNamespace("xcms", quietly = TRUE)) {
  stop("Package 'xcms' is required.")
}

if (!requireNamespace("CAMERA", quietly = TRUE)) {
  stop(
    paste(
      "Package 'CAMERA' is required for annotation.",
      "Please install CAMERA first."
    )
  )
}

# XCMS の中間結果を保存したフォルダ
input_dir <- "./output/xcms_pca"
output_dir <- "./output/xcms_annotation"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 事前に作成した xcms オブジェクトを読み込む
load(file.path(input_dir, "xcms_object.rds"))

if (!exists("xset")) {
  stop("xcms object 'xset' was not found in xcms_object.rds.")
}
```

### CAMERA xsAnnotate オブジェクトを作り、ピークをグループ化する

``` r
# 1. CAMERA xsAnnotate オブジェクトを作り、ピークをグループ化する
print("Creating CAMERA xsAnnotate object...")
an <- CAMERA::xsAnnotate(xset)

print("Grouping peaks with groupFWHM...")
an <- CAMERA::groupFWHM(an, perfwhm = 0.6)
```

### 同位体パターン候補を推定する

``` r
# 2. 同位体パターン候補を推定する
print("Finding isotope annotations...")
an <- CAMERA::findIsotopes(an, mzabs = 0.01, ppm = 5)
```

### 正イオンモードのアダクト候補を付ける

``` r
# 3. 正イオンモードのアダクト候補を付ける
print("Finding adduct annotations in positive mode...")
an <- CAMERA::findAdducts(an, polarity = "positive")
```

### annotation 結果の一覧を data.frame として取り出す

``` r
# 4. annotation 結果の一覧を data.frame として取り出す
annotation_table <- CAMERA::getPeaklist(an)
annotation_table <- as.data.frame(annotation_table)
```

### 結果を保存する

``` r
# 5. 結果を保存する
save(an, file = file.path(output_dir, "camera_annotation_object.rds"))
save(annotation_table, file = file.path(output_dir, "camera_annotation_table.rds"))
write.csv(
  annotation_table,
  file = file.path(output_dir, "camera_annotation_table.csv"),
  row.names = FALSE
)

print(paste("Annotated", nrow(annotation_table), "features"))
print(paste("Saved outputs to", output_dir))
```

## 3. 脂質候補として優先して見るピークを絞る

CAMERA の annotation をもとに、脂質解析に適した feature だけを残し、
次の MS1 照合に回す候補一覧を作ります。

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# このスクリプトでは、CAMERA の annotation 結果から
# リピドミクスで使いたい feature を前処理し、
# LipidBlast の full index と precursor m/z ベースで照合する。
#
# ここでは MS2 スペクトル照合までは行わず、
# まずは MS1 レベルの脂質候補一覧を作る。

annotation_file <- "./output/xcms_annotation/camera_annotation_table.csv"
library_index_file <- "./output/lipidblast_index/lipidblast_full_index.csv"
output_dir <- "./output/lipidblast_ms1_candidates"
mz_tolerance <- 0.01

if (!file.exists(annotation_file)) {
  stop("camera_annotation_table.csv was not found.")
}

if (!file.exists(library_index_file)) {
  stop("lipidblast_full_index.csv was not found.")
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

annotation_table <- read.csv(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
library_index <- read.csv(library_index_file, stringsAsFactors = FALSE)

annotation_table$feature_id <- paste0("feature_", seq_len(nrow(annotation_table)))
annotation_table$isotopes[is.na(annotation_table$isotopes)] <- ""
annotation_table$adduct[is.na(annotation_table$adduct)] <- ""
```

### モノアイソトープ feature とアダクト候補の判定関数を定義する

``` r
# 1. モノアイソトープ feature とアダクト候補の判定関数を定義する
is_monoisotope_like <- function(x) {
  if (is.na(x) || x == "") return(TRUE)
  grepl("\\[M\\](\\+|[0-9]\\+)$", x)
}

# リピドミクスでまず残したい主要アダクト
# 空欄は「annotation 未確定」なので候補集合には残す
is_major_lipid_adduct <- function(x) {
  if (is.na(x) || x == "") return(TRUE)
  grepl("^\\[M\\]\\+$", x) ||
    grepl("^\\[M\\+H\\]\\+", x) ||
    grepl("^\\[M\\+NH4\\]\\+", x) ||
    grepl("^\\[M\\+Na\\]\\+", x) ||
    grepl("^\\[M\\+K\\]\\+", x) ||
    grepl("^\\[M\\+2H\\]2\\+", x) ||
    grepl("^\\[M\\+3H\\]3\\+", x)
}
```

### フィルタリングを適用して脂質候補 feature を絞り込み、保存する

``` r
# 2. フィルタリングを適用して脂質候補 feature を絞り込み、保存する
annotation_table$is_monoisotope_like <- vapply(annotation_table$isotopes, is_monoisotope_like, logical(1))
annotation_table$is_major_lipid_adduct <- vapply(annotation_table$adduct, is_major_lipid_adduct, logical(1))

candidate_features <- annotation_table[
  annotation_table$is_monoisotope_like & annotation_table$is_major_lipid_adduct,
  ,
  drop = FALSE
]

write.csv(
  candidate_features,
  file = file.path(output_dir, "camera_ms1_candidate_features.csv"),
  row.names = FALSE
)
```

### LipidBlast インデックスを m/z でソートして feature と照合する

``` r
# 3. LipidBlast インデックスを m/z でソートして feature と照合する
ord <- order(library_index$PRECURSORMZ)
library_index <- library_index[ord, , drop = FALSE]
library_mz <- library_index$PRECURSORMZ

match_rows <- vector("list", nrow(candidate_features))

for (i in seq_len(nrow(candidate_features))) {
  target_mz <- candidate_features$mz[i]
  lower <- target_mz - mz_tolerance
  upper <- target_mz + mz_tolerance

  start_idx <- findInterval(lower, library_mz) + 1L
  end_idx <- findInterval(upper, library_mz)

  if (start_idx > end_idx || start_idx < 1L || end_idx < 1L) {
    next
  }

  hit <- library_index[start_idx:end_idx, , drop = FALSE]
  if (nrow(hit) == 0) {
    next
  }

  match_rows[[i]] <- data.frame(
    feature_id = candidate_features$feature_id[i],
    feature_mz = candidate_features$mz[i],
    feature_rt = candidate_features$rt[i],
    camera_adduct = candidate_features$adduct[i],
    camera_isotopes = candidate_features$isotopes[i],
    camera_pcgroup = candidate_features$pcgroup[i],
    mz_diff = candidate_features$mz[i] - hit$PRECURSORMZ,
    library_entry_id = hit$entry_id,
    library_name = hit$NAME,
    library_precursor_mz = hit$PRECURSORMZ,
    library_precursor_type = hit$PRECURSORTYPE,
    library_formula = hit$FORMULA,
    library_inchikey = hit$INCHIKEY,
    library_compound_class = hit$COMPOUNDCLASS,
    library_retention_time = hit$RETENTIONTIME,
    library_start_line = hit$start_line,
    library_end_line = hit$end_line,
    stringsAsFactors = FALSE
  )
}
```

### 照合結果とサマリーを保存する

``` r
# 4. 照合結果とサマリーを保存する
match_rows <- Filter(Negate(is.null), match_rows)
ms1_matches <- if (length(match_rows) == 0) data.frame() else do.call(rbind, match_rows)

write.csv(
  ms1_matches,
  file = file.path(output_dir, "camera_lipidblast_ms1_matches.csv"),
  row.names = FALSE
)
save(
  ms1_matches,
  file = file.path(output_dir, "camera_lipidblast_ms1_matches.rds")
)

summary_table <- data.frame(
  total_camera_features = nrow(annotation_table),
  filtered_camera_features = nrow(candidate_features),
  matched_features = if (nrow(ms1_matches) == 0) 0 else length(unique(ms1_matches$feature_id)),
  candidate_match_rows = nrow(ms1_matches),
  mz_tolerance = mz_tolerance,
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file = file.path(output_dir, "camera_lipidblast_ms1_matches_summary.csv"),
  row.names = FALSE
)

print(paste("Total CAMERA features:", nrow(annotation_table)))
print(paste("Filtered CAMERA features:", nrow(candidate_features)))
print(paste("Matched features:", if (nrow(ms1_matches) == 0) 0 else length(unique(ms1_matches$feature_id))))
print(paste("Candidate match rows:", nrow(ms1_matches)))
print(paste("Saved outputs to", output_dir))
```

## 4. MS1 の質量から脂質候補名を付ける

LIPID MAPS の構造データベースを検索用 DB に変換し、 各 feature の m/z
から一致する脂質候補名を探します。

### LIPID MAPS を検索用 DB に変換する

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# LIPID MAPS LMSD の structures.sdf から CompDb を作成する。

if (!requireNamespace("CompoundDb", quietly = TRUE)) {
  stop("Package 'CompoundDb' is required. Please install it before running this script.")
}

library(CompoundDb)

sdf_file <- "./output/lipidmaps/structures.sdf"
output_dir <- "./output/lipidmaps"
db_file_name <- "CompDb.LIPIDMAPS.LMSD.2026-03-18.sqlite"

if (!file.exists(sdf_file)) {
  stop("structures.sdf was not found.")
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
```

### LMSD SDF から化合物テーブルを読み込む

``` r
# 1. LMSD SDF から化合物テーブルを読み込む
message("Reading LMSD SDF ...")
cmps <- compound_tbl_sdf(sdf_file)
```

### メタデータを作り CompDb ファイルを生成する

``` r
# 2. メタデータを作り CompDb ファイルを生成する
message("Preparing metadata ...")
metad <- make_metadata(
  source = "LIPID MAPS LMSD",
  url = "https://www.lipidmaps.org/databases/lmsd",
  source_version = "2026-03-18",
  source_date = "2026-03-18",
  organism = NA_character_
)

message("Creating CompDb ...")
db_file <- createCompDb(
  cmps,
  metadata = metad,
  path = output_dir,
  dbFile = db_file_name
)

cdb <- CompDb(db_file)
```

### 作成した CompDb からサマリーと preview を保存する

``` r
# 3. 作成した CompDb からサマリーと preview を保存する
summary_table <- data.frame(
  db_file = db_file,
  compound_count = nrow(compounds(cdb, columns = "compound_id")),
  stringsAsFactors = FALSE
)

preview_cols <- intersect(
  c("compound_id", "name", "formula", "exactmass", "inchikey", "smiles"),
  colnames(compounds(cdb))
)

preview_table <- compounds(cdb, columns = preview_cols)
if (nrow(preview_table) > 20) {
  preview_table <- preview_table[seq_len(20), , drop = FALSE]
}

write.csv(
  summary_table,
  file = file.path(output_dir, "lipidmaps_compdb_summary.csv"),
  row.names = FALSE
)

write.csv(
  preview_table,
  file = file.path(output_dir, "lipidmaps_compdb_preview.csv"),
  row.names = FALSE
)

message("CompDb saved to: ", db_file)
message("Compound count: ", summary_table$compound_count)
```

### CompDb を使って m/z を照合し、MS1 ベースの候補名を付ける

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# LIPID MAPS CompDb を使って、MS1 feature に対する脂質候補を付ける。
# 速度を優先するため、CompDb から compound 情報だけを一度読み出し、
# 主要アダクトごとに理論 m/z を前計算して照合する。

if (!requireNamespace("CompoundDb", quietly = TRUE)) {
  stop("Package 'CompoundDb' is required.")
}

library(CompoundDb)

feature_file <- "./output/lipidblast_ms1_candidates/camera_ms1_candidate_features.csv"
compdb_file <- "./output/lipidmaps/CompDb.LIPIDMAPS.LMSD.2026-03-18.sqlite"
output_dir <- "./output/lipidmaps_ms1_candidates"

mz_tolerance <- 0.01
top_n_per_feature <- 10

if (!file.exists(feature_file)) stop("camera_ms1_candidate_features.csv was not found.")
if (!file.exists(compdb_file)) stop("CompDb file was not found.")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
```

### アダクト定義と feature 表を準備する

``` r
# 1. アダクト定義と feature 表を準備する
adduct_definitions <- data.frame(
  adduct = c("[M]+", "[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+", "[M+2H]2+", "[M+3H]3+"),
  delta_mass = c(0, 1.007276466812, 18.033823, 22.989218, 38.963158, 2 * 1.007276466812, 3 * 1.007276466812),
  charge = c(1, 1, 1, 1, 1, 2, 3),
  stringsAsFactors = FALSE
)

feature_table <- read.csv(feature_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"feature_id" %in% names(feature_table)) {
  stop("feature_id column was not found.")
}

# CAMERA の adduct 列が空なら、主要アダクトを広めに見る。
feature_table$ms1_query_adducts <- vapply(seq_len(nrow(feature_table)), function(i) {
  adduct <- trimws(feature_table$adduct[i])
  if (identical(adduct, "") || is.na(adduct)) {
    paste(adduct_definitions$adduct, collapse = "|")
  } else if (adduct %in% adduct_definitions$adduct) {
    adduct
  } else {
    paste(adduct_definitions$adduct, collapse = "|")
  }
}, FUN.VALUE = character(1))
```

### CompDb から化合物情報を読み込み、アダクトごとに理論 m/z を前計算する

``` r
# 2. CompDb から化合物情報を読み込み、アダクトごとに理論 m/z を前計算する
cdb <- CompDb(compdb_file)
compound_table <- compounds(
  cdb,
  columns = c("compound_id", "name", "formula", "exactmass", "inchikey", "smiles")
)

compound_table <- as.data.frame(compound_table, stringsAsFactors = FALSE)
compound_table <- compound_table[!is.na(compound_table$exactmass), , drop = FALSE]
compound_table$lipidmaps_category_code <- substr(compound_table$compound_id, 1, 4)

match_rows <- list()
match_index <- 0L
```

### アダクトごとに理論 m/z と feature m/z を照合する

``` r
# 3. アダクトごとに理論 m/z と feature m/z を照合する
for (k in seq_len(nrow(adduct_definitions))) {
  adduct_name <- adduct_definitions$adduct[k]
  delta_mass <- adduct_definitions$delta_mass[k]
  charge <- adduct_definitions$charge[k]

  target_df <- compound_table
  target_df$query_adduct <- adduct_name
  target_df$expected_mz <- (target_df$exactmass + delta_mass) / charge
  target_df <- target_df[order(target_df$expected_mz), , drop = FALSE]

  feature_idx <- grepl(adduct_name, feature_table$ms1_query_adducts, fixed = TRUE)
  query_df <- feature_table[feature_idx, c("feature_id", "mz", "rt", "adduct"), drop = FALSE]
  if (nrow(query_df) == 0) next

  for (i in seq_len(nrow(query_df))) {
    mz_value <- query_df$mz[i]
    lower <- mz_value - mz_tolerance
    upper <- mz_value + mz_tolerance

    left <- findInterval(lower, target_df$expected_mz)
    right <- findInterval(upper, target_df$expected_mz)
    if (right <= left) next

    hit_rows <- target_df[(left + 1):right, , drop = FALSE]
    if (nrow(hit_rows) == 0) next

    hit_rows$feature_id <- query_df$feature_id[i]
    hit_rows$feature_mz <- mz_value
    hit_rows$feature_rt <- query_df$rt[i]
    hit_rows$feature_camera_adduct <- query_df$adduct[i]
    hit_rows$mz_diff <- hit_rows$expected_mz - mz_value
    hit_rows$abs_mz_diff <- abs(hit_rows$mz_diff)
    hit_rows$ppm_error <- (hit_rows$mz_diff / mz_value) * 1e6
    hit_rows$abs_ppm_error <- abs(hit_rows$ppm_error)

    match_index <- match_index + 1L
    match_rows[[match_index]] <- hit_rows[, c(
      "feature_id",
      "feature_mz",
      "feature_rt",
      "feature_camera_adduct",
      "query_adduct",
      "mz_diff",
      "abs_mz_diff",
      "ppm_error",
      "abs_ppm_error",
      "compound_id",
      "name",
      "formula",
      "exactmass",
      "inchikey",
      "smiles",
      "lipidmaps_category_code"
    )]
  }
}
```

### 上位候補と feature ごとの最良候補を選ぶ

``` r
# 4. 上位候補と feature ごとの最良候補を選ぶ
all_matches <- if (length(match_rows) == 0) data.frame() else do.call(rbind, match_rows)

if (nrow(all_matches) > 0) {
  all_matches <- all_matches[order(all_matches$feature_id, all_matches$abs_ppm_error, all_matches$query_adduct, all_matches$compound_id), , drop = FALSE]
  best_matches <- do.call(
    rbind,
    lapply(split(all_matches, all_matches$feature_id), function(df) head(df, top_n_per_feature))
  )
} else {
  best_matches <- data.frame()
}

feature_best_match <- if (nrow(best_matches) == 0) {
  data.frame()
} else {
  do.call(
    rbind,
    lapply(split(best_matches, best_matches$feature_id), function(df) df[1, , drop = FALSE])
  )
}
```

### 照合結果をファイルに保存する

``` r
# 5. 照合結果をファイルに保存する
summary_table <- data.frame(
  total_candidate_features = nrow(feature_table),
  matched_feature_count = if (nrow(all_matches) == 0) 0 else length(unique(all_matches$feature_id)),
  total_match_rows = nrow(all_matches),
  top_match_rows = nrow(best_matches),
  mz_tolerance = mz_tolerance,
  ranking_metric = "abs_ppm_error",
  stringsAsFactors = FALSE
)

write.csv(
  all_matches,
  file = file.path(output_dir, "lipidmaps_ms1_all_matches.csv"),
  row.names = FALSE
)

write.csv(
  best_matches,
  file = file.path(output_dir, "lipidmaps_ms1_top_matches.csv"),
  row.names = FALSE
)

write.csv(
  feature_best_match,
  file = file.path(output_dir, "lipidmaps_ms1_feature_best_matches.csv"),
  row.names = FALSE
)

write.csv(
  summary_table,
  file = file.path(output_dir, "lipidmaps_ms1_match_summary.csv"),
  row.names = FALSE
)

print(paste("Matched features:", summary_table$matched_feature_count))
print(paste("Total match rows:", summary_table$total_match_rows))
print(paste("Top match rows:", summary_table$top_match_rows))
```

## 5. MS2 ライブラリ照合の準備をする

MS1 の上位候補をもとに LipidBlast から必要なスペクトルだけを抽出し、 MS2
照合に使うライブラリを絞り込みます。

### MS1 上位候補から LipidBlast 用の request table を作る

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# MS1 候補表から、feature ごとに mz 差が小さい上位候補だけを残す。

input_file <- "./output/lipidblast_ms1_candidates/camera_lipidblast_ms1_matches.csv"
output_dir <- "./output/lipidblast_ms2_candidates"
top_n <- 5

if (!file.exists(input_file)) {
  stop("camera_lipidblast_ms1_matches.csv was not found.")
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
```

### feature ごとに m/z 差が小さい上位候補を絞り込む

``` r
# 1. feature ごとに m/z 差が小さい上位候補を絞り込む
x <- read.csv(input_file, stringsAsFactors = FALSE)
x$abs_mz_diff <- abs(x$mz_diff)
x <- x[order(x$feature_id, x$abs_mz_diff, x$library_entry_id), , drop = FALSE]

split_list <- split(x, x$feature_id)
top_candidates <- do.call(
  rbind,
  lapply(split_list, function(df) head(df, top_n))
)

top_candidates <- top_candidates[order(top_candidates$feature_id, top_candidates$abs_mz_diff), , drop = FALSE]

entry_ids <- data.frame(
  library_entry_id = sort(unique(top_candidates$library_entry_id)),
  stringsAsFactors = FALSE
)

summary_table <- data.frame(
  top_n = top_n,
  matched_feature_count = length(unique(top_candidates$feature_id)),
  candidate_row_count = nrow(top_candidates),
  unique_library_entry_count = nrow(entry_ids),
  stringsAsFactors = FALSE
)
```

### 結果をファイルに保存する

``` r
# 2. 結果をファイルに保存する
write.csv(
  top_candidates,
  file = file.path(output_dir, "feature_top_lipidblast_candidates.csv"),
  row.names = FALSE
)

write.csv(
  entry_ids,
  file = file.path(output_dir, "top_lipidblast_entry_ids.csv"),
  row.names = FALSE
)

write.csv(
  summary_table,
  file = file.path(output_dir, "feature_top_lipidblast_candidates_summary.csv"),
  row.names = FALSE
)

print(paste("Matched features:", summary_table$matched_feature_count))
print(paste("Candidate rows:", summary_table$candidate_row_count))
print(paste("Unique library entries:", summary_table$unique_library_entry_count))
```

### request table の候補だけを LipidBlast から抽出する

### 入力ファイル・出力先・基本設定

``` powershell
param(
  [string]$MspFile = "C:\Users\yamamoto\R\MSplusR\dev\MSDIAL-TandemMassSpectralAtlas-VS69-Pos.msp",
  [string]$EntryIdCsv = "C:\Users\yamamoto\R\MSplusR\dev\lipidblast_ms2_candidates\top_lipidblast_entry_ids.csv",
  [string]$OutputDir = "C:\Users\yamamoto\R\MSplusR\dev\lipidblast_ms2_candidates"
)

if (-not (Test-Path $MspFile)) {
  throw "MSP file was not found: $MspFile"
}

if (-not (Test-Path $EntryIdCsv)) {
  throw "Entry ID CSV was not found: $EntryIdCsv"
}

New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null

$entryIdTable = Import-Csv $EntryIdCsv
$wanted = New-Object 'System.Collections.Generic.HashSet[int]'
foreach ($row in $entryIdTable) {
  [void]$wanted.Add([int]$row.library_entry_id)
}

$entries = New-Object System.Collections.Generic.List[object]
$peaks = New-Object System.Collections.Generic.List[object]

$reader = [System.IO.StreamReader]::new($MspFile)

$currentEntryId = 0
$capture = $false
$current = $null

function Save-CurrentEntry {
  param(
    [ref]$Current,
    [System.Collections.Generic.List[object]]$Entries,
    [System.Collections.Generic.List[object]]$Peaks
  )

  if ($null -eq $Current.Value) { return }

  $Entries.Add([pscustomobject]@{
    entry_id = $Current.Value.entry_id
    NAME = $Current.Value.NAME
    PRECURSORMZ = $Current.Value.PRECURSORMZ
    PRECURSORTYPE = $Current.Value.PRECURSORTYPE
    FORMULA = $Current.Value.FORMULA
    INCHIKEY = $Current.Value.INCHIKEY
    COMPOUNDCLASS = $Current.Value.COMPOUNDCLASS
    RETENTIONTIME = $Current.Value.RETENTIONTIME
    NUM_PEAKS = $Current.Value.NUM_PEAKS
    extracted_peak_count = $Current.Value.extracted_peak_count
  }) | Out-Null

  foreach ($peak in $Current.Value.peaks) {
    $Peaks.Add([pscustomobject]@{
      entry_id = $Current.Value.entry_id
      peak_mz = $peak.peak_mz
      peak_intensity = $peak.peak_intensity
    }) | Out-Null
  }
}

try {
  while (($line = $reader.ReadLine()) -ne $null) {
    if ($line.StartsWith("NAME:")) {
      if ($capture) {
        Save-CurrentEntry -Current ([ref]$current) -Entries $entries -Peaks $peaks
      }

      $currentEntryId += 1
      $capture = $wanted.Contains($currentEntryId)

      if ($capture) {
        $current = @{
          entry_id = $currentEntryId
          NAME = ($line.Substring(5)).Trim()
          PRECURSORMZ = $null
          PRECURSORTYPE = ""
          FORMULA = ""
          INCHIKEY = ""
          COMPOUNDCLASS = ""
          RETENTIONTIME = $null
          NUM_PEAKS = $null
          extracted_peak_count = 0
          peaks = New-Object System.Collections.Generic.List[object]
        }
      } else {
        $current = $null
      }
      continue
    }

    if (-not $capture) { continue }

    if ([string]::IsNullOrWhiteSpace($line)) { continue }

    if ($line -match '^PRECURSORMZ:\s*(.+)$') {
      $current.PRECURSORMZ = [double]$matches[1]
      continue
    }
    if ($line -match '^PRECURSORTYPE:\s*(.+)$') {
      $current.PRECURSORTYPE = $matches[1].Trim()
      continue
    }
    if ($line -match '^FORMULA:\s*(.+)$') {
      $current.FORMULA = $matches[1].Trim()
      continue
    }
    if ($line -match '^INCHIKEY:\s*(.+)$') {
      $current.INCHIKEY = $matches[1].Trim()
      continue
    }
    if ($line -match '^COMPOUNDCLASS:\s*(.+)$') {
      $current.COMPOUNDCLASS = $matches[1].Trim()
      continue
    }
    if ($line -match '^RETENTIONTIME:\s*(.+)$') {
      $current.RETENTIONTIME = [double]$matches[1]
      continue
    }
    if ($line -match '^Num Peaks:\s*(.+)$') {
      $current.NUM_PEAKS = [int]$matches[1]
      continue
    }
    if ($line -match '^\s*([0-9]+(?:\.[0-9]+)?)\s+([0-9]+(?:\.[0-9]+)?)\s*$') {
      $current.peaks.Add([pscustomobject]@{
        peak_mz = [double]$matches[1]
        peak_intensity = [double]$matches[2]
      }) | Out-Null
      $current.extracted_peak_count += 1
    }
  }

  if ($capture) {
    Save-CurrentEntry -Current ([ref]$current) -Entries $entries -Peaks $peaks
  }
}
finally {
  $reader.Close()
}

$entryOut = Join-Path $OutputDir "top_lipidblast_entry_table.csv"
$peakOut = Join-Path $OutputDir "top_lipidblast_peak_table.csv"
$summaryOut = Join-Path $OutputDir "top_lipidblast_extract_summary.csv"

$entries | Export-Csv -Path $entryOut -NoTypeInformation -Encoding UTF8
$peaks | Export-Csv -Path $peakOut -NoTypeInformation -Encoding UTF8

[pscustomobject]@{
  requested_entry_count = $wanted.Count
  extracted_entry_count = $entries.Count
  extracted_peak_row_count = $peaks.Count
} | Export-Csv -Path $summaryOut -NoTypeInformation -Encoding UTF8

Write-Output "Requested entries: $($wanted.Count)"
Write-Output "Extracted entries: $($entries.Count)"
Write-Output "Extracted peak rows: $($peaks.Count)"
```

## 6. 実測 MS2 をピーク候補ごとにまとめる

実サンプルの mzML から DDA MS2 を取り出し、 feature ごとに consensus
spectrum を作ります。

### 実測 DDA MS2 を feature ごとに整理する

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# このスクリプトでは、MS1 で LipidBlast 候補が付いた feature に対して、
# 実サンプルの mzML から DDA MS2 スペクトルを抽出する。
#
# feature と MS2 スペクトルの紐付けは、以下の条件で行う。
# - precursor m/z が feature m/z に近い
# - retention time が feature の rtmin-rtmax 範囲に入る

if (!requireNamespace("MSnbase", quietly = TRUE)) {
  stop("Package 'MSnbase' is required.")
}

library(MSnbase)

feature_file <- "./output/lipidblast_ms1_candidates/camera_ms1_candidate_features.csv"
ms1_match_file <- "./output/lipidblast_ms1_candidates/camera_lipidblast_ms1_matches.csv"
mzml_dir <- "./mzml"
output_dir <- "./output/dda_feature_spectra"

precursor_mz_tolerance <- 0.01
rt_padding <- 5

if (!file.exists(feature_file)) {
  stop("camera_ms1_candidate_features.csv was not found.")
}

if (!file.exists(ms1_match_file)) {
  stop("camera_lipidblast_ms1_matches.csv was not found.")
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

feature_table <- read.csv(feature_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!"feature_id" %in% names(feature_table)) {
  feature_table$feature_id <- paste0("feature_", seq_len(nrow(feature_table)))
}

ms1_matches <- read.csv(ms1_match_file, stringsAsFactors = FALSE)
matched_feature_ids <- unique(ms1_matches$feature_id)

feature_table <- feature_table[feature_table$feature_id %in% matched_feature_ids, , drop = FALSE]

mzfiles <- list.files(mzml_dir, pattern = "\\.mzML$", full.names = TRUE)
if (length(mzfiles) == 0) {
  stop("No mzML files were found.")
}

match_rows <- list()
matched_spectra <- list()
match_index <- 0L
```

### mzML を順に読み込んで DDA MS2 スペクトルを feature と照合する

``` r
# 1. mzML を順に読み込んで DDA MS2 スペクトルを feature と照合する
for (mzfile in mzfiles) {
  sample_name <- tools::file_path_sans_ext(basename(mzfile))
  print(paste("Processing sample:", sample_name))

  if (sample_name %in% names(feature_table)) {
    sample_feature_table <- feature_table[feature_table[[sample_name]] > 0, , drop = FALSE]
  } else {
    sample_feature_table <- feature_table
  }

  if (nrow(sample_feature_table) == 0) next

  x <- readMSData(mzfile, mode = "onDisk")
  h <- fData(x)

  idx2 <- which(h$msLevel == 2 & !is.na(h$precursorMZ))
  if (length(idx2) == 0) next

  ms2 <- x[idx2]
  h2 <- h[idx2, , drop = FALSE]

  spectra_list <- spectra(ms2)

  for (i in seq_len(nrow(sample_feature_table))) {
    feature_id <- sample_feature_table$feature_id[i]
    feature_mz <- sample_feature_table$mz[i]
    lower_mz <- feature_mz - precursor_mz_tolerance
    upper_mz <- feature_mz + precursor_mz_tolerance
    lower_rt <- sample_feature_table$rtmin[i] - rt_padding
    upper_rt <- sample_feature_table$rtmax[i] + rt_padding

    hit_idx <- which(
      h2$precursorMZ >= lower_mz &
        h2$precursorMZ <= upper_mz &
        h2$retentionTime >= lower_rt &
        h2$retentionTime <= upper_rt
    )

    if (length(hit_idx) == 0) next

    for (j in hit_idx) {
      sp <- spectra_list[[j]]
      if (length(sp@mz) == 0) next

      match_index <- match_index + 1L

        match_rows[[match_index]] <- data.frame(
        match_id = match_index,
        feature_id = feature_id,
        sample_name = sample_name,
        feature_mz = feature_mz,
        feature_rt = sample_feature_table$rt[i],
        feature_rtmin = sample_feature_table$rtmin[i],
        feature_rtmax = sample_feature_table$rtmax[i],
        spectrum_precursor_mz = h2$precursorMZ[j],
        spectrum_rt = h2$retentionTime[j],
        acquisition_num = h2$acquisitionNum[j],
        peaks_count = length(sp@mz),
        stringsAsFactors = FALSE
      )

      matched_spectra[[match_index]] <- list(
        match_id = match_index,
        feature_id = feature_id,
        sample_name = sample_name,
        precursor_mz = h2$precursorMZ[j],
        retention_time = h2$retentionTime[j],
        mz = sp@mz,
        intensity = sp@intensity
      )
    }
  }
}
```

### 照合結果とスペクトルをファイルに保存する

``` r
# 2. 照合結果とスペクトルをファイルに保存する
dda_match_table <- if (length(match_rows) == 0) data.frame() else do.call(rbind, match_rows)

write.csv(
  dda_match_table,
  file = file.path(output_dir, "dda_feature_match_table.csv"),
  row.names = FALSE
)
save(
  dda_match_table,
  file = file.path(output_dir, "dda_feature_match_table.rds")
)
save(
  matched_spectra,
  file = file.path(output_dir, "dda_feature_spectra.rds")
)

summary_table <- data.frame(
  matched_feature_count = if (nrow(dda_match_table) == 0) 0 else length(unique(dda_match_table$feature_id)),
  matched_spectrum_count = nrow(dda_match_table),
  sample_count = if (nrow(dda_match_table) == 0) 0 else length(unique(dda_match_table$sample_name)),
  precursor_mz_tolerance = precursor_mz_tolerance,
  rt_padding = rt_padding,
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file = file.path(output_dir, "dda_feature_match_summary.csv"),
  row.names = FALSE
)

print(paste("Matched features:", summary_table$matched_feature_count))
print(paste("Matched spectra:", summary_table$matched_spectrum_count))
print(paste("Saved outputs to", output_dir))
```

### 複数 MS2 を consensus spectrum にまとめる

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# DDA から抽出した複数の MS2 スペクトルを、
# feature_id x sample_name ごとに 1 本の consensus spectrum に統合する。

if (!requireNamespace("MSnbase", quietly = TRUE)) {
  stop("Package 'MSnbase' is required.")
}

library(MSnbase)

input_match_table <- "./output/dda_feature_spectra/dda_feature_match_table.csv"
input_spectra_rds <- "./output/dda_feature_spectra/dda_feature_spectra.rds"
output_dir <- "./output/dda_feature_spectra"

mzd_tolerance <- 0.01
min_peak_fraction <- 0.3

if (!file.exists(input_match_table)) {
  stop("dda_feature_match_table.csv was not found.")
}

if (!file.exists(input_spectra_rds)) {
  stop("dda_feature_spectra.rds was not found.")
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

match_table <- read.csv(input_match_table, stringsAsFactors = FALSE)
load(input_spectra_rds)

if (!exists("matched_spectra")) {
  stop("Object 'matched_spectra' was not found in dda_feature_spectra.rds")
}

if (length(matched_spectra) == 0 || nrow(match_table) == 0) {
  stop("No matched spectra were available for consensus generation.")
}

group_id <- paste(match_table$feature_id, match_table$sample_name, sep = "||")
match_table$group_id <- group_id
```

### 各スペクトルを Spectrum2 オブジェクトに変換する

``` r
# 1. 各スペクトルを Spectrum2 オブジェクトに変換する
spectrum_objects <- lapply(seq_along(matched_spectra), function(i) {
  sp <- matched_spectra[[i]]
  MSnbase:::Spectrum2(
    mz = as.double(sp$mz),
    intensity = as.double(sp$intensity),
    peaksCount = as.integer(length(sp$mz)),
    rt = as.double(sp$retention_time),
    acquisitionNum = as.integer(sp$match_id),
    precursorMz = as.double(sp$precursor_mz)
  )
})

element_metadata <- match_table[, c(
  "match_id",
  "feature_id",
  "sample_name",
  "feature_mz",
  "feature_rt",
  "feature_rtmin",
  "feature_rtmax",
  "spectrum_precursor_mz",
  "spectrum_rt",
  "acquisition_num",
  "peaks_count",
  "group_id"
)]
```

### MSpectra オブジェクトを組み立てて consensus spectrum を計算する

``` r
# 2. MSpectra オブジェクトを組み立てて consensus spectrum を計算する
spectra_obj <- MSpectra(
  spectrum_objects,
  elementMetadata = S4Vectors::DataFrame(element_metadata)
)

consensus_obj <- combineSpectra(
  spectra_obj,
  fcol = "group_id",
  method = consensusSpectrum,
  mzd = mzd_tolerance,
  minProp = min_peak_fraction,
  intensityFun = stats::median,
  mzFun = stats::median
)

consensus_meta <- as.data.frame(S4Vectors::mcols(consensus_obj), stringsAsFactors = FALSE)
consensus_meta$consensus_peaks_count <- peaksCount(consensus_obj)
consensus_meta$consensus_precursor_mz <- precursorMz(consensus_obj)
consensus_meta$consensus_rt <- rtime(consensus_obj)
consensus_meta$source_spectra_count <- as.integer(table(group_id)[consensus_meta$group_id])
```

### consensus スペクトルのメタ情報をリストに整理する

``` r
# 3. consensus スペクトルのメタ情報をリストに整理する
consensus_spectra <- lapply(seq_along(consensus_obj), function(i) {
  list(
    group_id = consensus_meta$group_id[i],
    feature_id = consensus_meta$feature_id[i],
    sample_name = consensus_meta$sample_name[i],
    feature_mz = consensus_meta$feature_mz[i],
    feature_rt = consensus_meta$feature_rt[i],
    source_spectra_count = consensus_meta$source_spectra_count[i],
    consensus_precursor_mz = consensus_meta$consensus_precursor_mz[i],
    consensus_rt = consensus_meta$consensus_rt[i],
    mz = mz(consensus_obj[[i]]),
    intensity = intensity(consensus_obj[[i]])
  )
})
```

### 結果をファイルに保存する

``` r
# 4. 結果をファイルに保存する
summary_table <- data.frame(
  consensus_group_count = length(consensus_spectra),
  original_spectrum_count = nrow(match_table),
  feature_count = length(unique(consensus_meta$feature_id)),
  sample_count = length(unique(consensus_meta$sample_name)),
  mzd_tolerance = mzd_tolerance,
  min_peak_fraction = min_peak_fraction,
  stringsAsFactors = FALSE
)

write.csv(
  consensus_meta,
  file = file.path(output_dir, "dda_consensus_match_table.csv"),
  row.names = FALSE
)

save(
  consensus_spectra,
  file = file.path(output_dir, "dda_consensus_spectra.rds")
)

save(
  consensus_obj,
  file = file.path(output_dir, "dda_consensus_mspectra.rds")
)

write.csv(
  summary_table,
  file = file.path(output_dir, "dda_consensus_summary.csv"),
  row.names = FALSE
)

print(paste("Consensus groups:", summary_table$consensus_group_count))
print(paste("Original spectra:", summary_table$original_spectrum_count))
print(paste("Features:", summary_table$feature_count))
print(paste("Saved outputs to", output_dir))
```

## 7. 実測 MS2 をライブラリと比べる

compareSpectra を使って consensus MS2 と LipidBlast スペクトルを比較し、
フラグメントパターンに基づく脂質候補を絞り込みます。

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# Spectra::compareSpectra を使って、consensus MS2 と
# LipidBlast 候補スペクトルを照合する。

if (!requireNamespace("Spectra", quietly = TRUE)) {
  stop("Package 'Spectra' is required.")
}

library(Spectra)

candidate_file <- "./output/lipidblast_ms2_candidates/feature_top_lipidblast_candidates.csv"
entry_file <- "./output/lipidblast_ms2_candidates/top_lipidblast_entry_table.csv"
peak_file <- "./output/lipidblast_ms2_candidates/top_lipidblast_peak_table.csv"
consensus_file <- "./output/dda_feature_spectra/dda_consensus_spectra.rds"
output_dir <- "./output/lipidblast_ms2_matching_comparespectra"

fragment_tolerance <- 0.02
min_query_peaks <- 3

if (!file.exists(candidate_file)) stop("feature_top_lipidblast_candidates.csv was not found.")
if (!file.exists(entry_file)) stop("top_lipidblast_entry_table.csv was not found.")
if (!file.exists(peak_file)) stop("top_lipidblast_peak_table.csv was not found.")
if (!file.exists(consensus_file)) stop("dda_consensus_spectra.rds was not found.")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

candidates <- read.csv(candidate_file, stringsAsFactors = FALSE)
entry_table <- read.csv(entry_file, stringsAsFactors = FALSE)
peak_table <- read.csv(peak_file, stringsAsFactors = FALSE)
load(consensus_file)

if (!exists("consensus_spectra")) {
  stop("Object 'consensus_spectra' was not found.")
}

consensus_spectra <- consensus_spectra[vapply(consensus_spectra, function(x) length(x$mz) >= min_query_peaks, logical(1))]

peak_split <- split(peak_table, peak_table$entry_id)
entry_table <- entry_table[match(unique(entry_table$entry_id), entry_table$entry_id), , drop = FALSE]
```

### LipidBlast 候補スペクトルの Spectra オブジェクトを作る

``` r
# 1. LipidBlast 候補スペクトルの Spectra オブジェクトを作る
target_df <- S4Vectors::DataFrame(
  msLevel = rep(2L, nrow(entry_table)),
  precursorMz = entry_table$PRECURSORMZ,
  rtime = entry_table$RETENTIONTIME,
  mz = I(lapply(as.character(entry_table$entry_id), function(id) peak_split[[id]]$peak_mz)),
  intensity = I(lapply(as.character(entry_table$entry_id), function(id) peak_split[[id]]$peak_intensity)),
  library_entry_id = entry_table$entry_id,
  library_name = entry_table$NAME,
  library_precursor_type = entry_table$PRECURSORTYPE,
  library_compound_class = entry_table$COMPOUNDCLASS
)
target_spectra <- Spectra(target_df, source = MsBackendDataFrame())
```

### consensus MS2 スペクトルの Spectra オブジェクトを作る

``` r
# 2. consensus MS2 スペクトルの Spectra オブジェクトを作る
query_df <- S4Vectors::DataFrame(
  msLevel = rep(2L, length(consensus_spectra)),
  precursorMz = vapply(consensus_spectra, function(x) x$consensus_precursor_mz, numeric(1)),
  rtime = vapply(consensus_spectra, function(x) x$consensus_rt, numeric(1)),
  mz = I(lapply(consensus_spectra, function(x) x$mz)),
  intensity = I(lapply(consensus_spectra, function(x) x$intensity)),
  group_id = vapply(consensus_spectra, function(x) x$group_id, character(1)),
  feature_id = vapply(consensus_spectra, function(x) x$feature_id, character(1)),
  sample_name = vapply(consensus_spectra, function(x) x$sample_name, character(1)),
  feature_mz = vapply(consensus_spectra, function(x) x$feature_mz, numeric(1)),
  feature_rt = vapply(consensus_spectra, function(x) x$feature_rt, numeric(1)),
  source_spectra_count = vapply(consensus_spectra, function(x) x$source_spectra_count, numeric(1))
)
query_spectra <- Spectra(query_df, source = MsBackendDataFrame())
```

### feature ごとに compareSpectra でスコアを計算する

``` r
# 3. feature ごとに compareSpectra でスコアを計算する
candidate_split <- split(candidates, candidates$feature_id)
query_indices <- split(seq_along(consensus_spectra), query_df$feature_id)
target_indices <- split(seq_len(nrow(entry_table)), entry_table$entry_id)

match_peak_count <- function(query_mz, target_mz, tolerance = 0.02) {
  if (length(query_mz) == 0 || length(target_mz) == 0) return(0L)
  used <- rep(FALSE, length(target_mz))
  count <- 0L
  for (i in seq_along(query_mz)) {
    diff <- abs(target_mz - query_mz[i])
    diff[used] <- Inf
    best <- which.min(diff)
    if (length(best) == 0 || is.infinite(diff[best]) || diff[best] > tolerance) next
    used[best] <- TRUE
    count <- count + 1L
  }
  count
}

result_rows <- list()
result_index <- 0L

for (feature_id in intersect(names(candidate_split), names(query_indices))) {
  feature_candidates <- candidate_split[[feature_id]]
  q_idx <- query_indices[[feature_id]]
  t_idx <- unlist(target_indices[as.character(feature_candidates$library_entry_id)], use.names = FALSE)
  t_idx <- unique(t_idx[!is.na(t_idx)])
  if (length(q_idx) == 0 || length(t_idx) == 0) next

  q_sp <- query_spectra[q_idx]
  t_sp <- target_spectra[t_idx]

  score_mat <- compareSpectra(
    q_sp,
    t_sp,
    tolerance = fragment_tolerance,
    ppm = 0
  )

  if (is.null(dim(score_mat))) {
    score_mat <- matrix(score_mat, nrow = length(q_sp), ncol = length(t_sp))
  }

  for (i in seq_len(nrow(score_mat))) {
    for (j in seq_len(ncol(score_mat))) {
      entry_id <- t_sp$library_entry_id[j]
      candidate_row <- feature_candidates[feature_candidates$library_entry_id == entry_id, , drop = FALSE]
      if (nrow(candidate_row) == 0) next

      matched_peak_n <- match_peak_count(
        mz(q_sp)[[i]],
        mz(t_sp)[[j]],
        tolerance = fragment_tolerance
      )

      result_index <- result_index + 1L
      result_rows[[result_index]] <- data.frame(
        group_id = q_sp$group_id[i],
        feature_id = q_sp$feature_id[i],
        sample_name = q_sp$sample_name[i],
        feature_mz = q_sp$feature_mz[i],
        feature_rt = q_sp$feature_rt[i],
        source_spectra_count = q_sp$source_spectra_count[i],
        consensus_peak_count = lengths(mz(q_sp))[i],
        library_entry_id = entry_id,
        library_name = t_sp$library_name[j],
        library_precursor_type = t_sp$library_precursor_type[j],
        library_compound_class = t_sp$library_compound_class[j],
        ms1_mz_diff = candidate_row$mz_diff[1],
        score = score_mat[i, j],
        matched_peak_count = matched_peak_n,
        stringsAsFactors = FALSE
      )
    }
  }
}
```

### 照合結果を整理してファイルに保存する

``` r
# 4. 照合結果を整理してファイルに保存する
match_table <- if (length(result_rows) == 0) data.frame() else do.call(rbind, result_rows)

if (nrow(match_table) > 0) {
  match_table <- match_table[order(
    match_table$group_id,
    -match_table$score,
    -match_table$matched_peak_count,
    abs(match_table$ms1_mz_diff)
  ), , drop = FALSE]
}

best_matches <- if (nrow(match_table) == 0) data.frame() else do.call(
  rbind,
  lapply(split(match_table, match_table$group_id), function(df) df[1, , drop = FALSE])
)

summary_table <- data.frame(
  consensus_group_count = length(consensus_spectra),
  compared_group_count = if (nrow(match_table) == 0) 0 else length(unique(match_table$group_id)),
  total_comparison_rows = nrow(match_table),
  best_match_rows = nrow(best_matches),
  fragment_tolerance = fragment_tolerance,
  min_query_peaks = min_query_peaks,
  stringsAsFactors = FALSE
)

write.csv(
  match_table,
  file = file.path(output_dir, "consensus_lipidblast_match_table_comparespectra.csv"),
  row.names = FALSE
)

write.csv(
  best_matches,
  file = file.path(output_dir, "consensus_lipidblast_best_matches_comparespectra.csv"),
  row.names = FALSE
)

write.csv(
  summary_table,
  file = file.path(output_dir, "consensus_lipidblast_match_summary_comparespectra.csv"),
  row.names = FALSE
)

print(paste("Compared groups:", summary_table$compared_group_count))
print(paste("Total comparison rows:", summary_table$total_comparison_rows))
print(paste("Best matches:", summary_table$best_match_rows))
```

## 8. MS1 と MS2 の結果を 1 つの annotation 表にまとめる

MS1 と MS2 それぞれの照合結果を結合し、 feature ごとの最終 annotation
matrix を作ります。

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

# MS1 は LIPID MAPS 候補を 1 本に絞らず保持し、
# MS2 は compareSpectra で得た代表候補を feature ごとに結合する。

feature_file <- "./output/lipidblast_ms1_candidates/camera_ms1_candidate_features.csv"
ms1_file <- "./output/lipidmaps_ms1_candidates/lipidmaps_ms1_all_matches.csv"
ms2_file <- "./output/lipidblast_ms2_matching_comparespectra/consensus_lipidblast_best_matches_comparespectra.csv"
output_dir <- "./output/lipid_annotation_matrix_combined_comparespectra"

collapse_unique <- function(x, limit = Inf) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(NA_character_)
  if (is.finite(limit)) x <- head(x, limit)
  paste(x, collapse = " | ")
}

if (!file.exists(feature_file)) stop("camera_ms1_candidate_features.csv was not found.")
if (!file.exists(ms1_file)) stop("lipidmaps_ms1_all_matches.csv was not found.")
if (!file.exists(ms2_file)) stop("consensus_lipidblast_best_matches_comparespectra.csv was not found.")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

feature_table <- read.csv(feature_file, stringsAsFactors = FALSE, check.names = FALSE)
ms1_all <- read.csv(ms1_file, stringsAsFactors = FALSE, check.names = FALSE)
ms2_best <- read.csv(ms2_file, stringsAsFactors = FALSE, check.names = FALSE)
```

### MS1 候補を feature ごとにまとめる

``` r
# 1. MS1 候補を feature ごとにまとめる
ms1_feature_summary <- do.call(
  rbind,
  lapply(split(ms1_all, ms1_all$feature_id), function(df) {
    df <- df[order(df$abs_ppm_error, df$query_adduct, df$compound_id), , drop = FALSE]
    data.frame(
      feature_id = df$feature_id[1],
      ms1_candidate_count = nrow(df),
      ms1_unique_compound_count = length(unique(df$compound_id)),
      ms1_query_adducts = collapse_unique(df$query_adduct),
      ms1_candidate_compound_ids = collapse_unique(df$compound_id, limit = 20),
      ms1_candidate_names = collapse_unique(df$name, limit = 20),
      ms1_candidate_formulas = collapse_unique(df$formula, limit = 20),
      ms1_category_codes = collapse_unique(df$lipidmaps_category_code),
      ms1_min_abs_mz_diff = min(df$abs_mz_diff, na.rm = TRUE),
      ms1_min_abs_ppm_error = min(df$abs_ppm_error, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
```

### MS2 候補を feature ごとにまとめて最良候補を選ぶ

``` r
# 2. MS2 候補を feature ごとにまとめて最良候補を選ぶ
ms2_feature_summary <- aggregate(
  cbind(
    ms2_sample_support = rep(1, nrow(ms2_best)),
    ms2_mean_score = ms2_best$score,
    ms2_max_score = ms2_best$score,
    ms2_mean_matched_peak_count = ms2_best$matched_peak_count,
    ms2_max_matched_peak_count = ms2_best$matched_peak_count
  ),
  by = list(
    feature_id = ms2_best$feature_id,
    ms2_library_entry_id = ms2_best$library_entry_id,
    ms2_annotation_name = ms2_best$library_name,
    ms2_annotation_precursor_type = ms2_best$library_precursor_type,
    ms2_annotation_compound_class = ms2_best$library_compound_class
  ),
  FUN = mean
)

ms2_feature_summary <- ms2_feature_summary[order(
  ms2_feature_summary$feature_id,
  -ms2_feature_summary$ms2_sample_support,
  -ms2_feature_summary$ms2_max_matched_peak_count,
  -ms2_feature_summary$ms2_mean_score,
  -ms2_feature_summary$ms2_max_score,
  ms2_feature_summary$ms2_library_entry_id
), , drop = FALSE]

ms2_feature_best <- do.call(
  rbind,
  lapply(split(ms2_feature_summary, ms2_feature_summary$feature_id), function(df) df[1, , drop = FALSE])
)
```

### MS1 と MS2 の結果を feature 表に結合して annotation 情報を付ける

``` r
# 3. MS1 と MS2 の結果を feature 表に結合して annotation 情報を付ける
combined <- merge(
  feature_table,
  ms1_feature_summary,
  by = "feature_id",
  all.x = TRUE,
  sort = FALSE
)
combined <- merge(combined, ms2_feature_best, by = "feature_id", all.x = TRUE, sort = FALSE)

combined$annotation_support_type <- ifelse(
  !is.na(combined$ms1_candidate_count) & !is.na(combined$ms2_library_entry_id),
  "MS1+MS2",
  ifelse(
    !is.na(combined$ms2_library_entry_id),
    "MS2_only",
    ifelse(!is.na(combined$ms1_candidate_count), "MS1_only", "none")
  )
)

combined$final_annotation_name <- ifelse(
  !is.na(combined$ms2_annotation_name),
  combined$ms2_annotation_name,
  NA
)

combined$final_annotation_source <- ifelse(
  !is.na(combined$ms2_annotation_name),
  "LipidBlast_MS2_compareSpectra",
  NA
)
```

### サマリー表と最終 annotation matrix をファイルに保存する

``` r
# 4. サマリー表と最終 annotation matrix をファイルに保存する
summary_table <- data.frame(
  total_candidate_features = nrow(combined),
  ms1_annotated_features = sum(!is.na(combined$ms1_candidate_count)),
  ms2_annotated_features = sum(!is.na(combined$ms2_library_entry_id)),
  ms1_ms2_annotated_features = sum(combined$annotation_support_type == "MS1+MS2"),
  ms1_annotation_mode = "retain_all_candidates",
  stringsAsFactors = FALSE
)

write.csv(
  combined,
  file = file.path(output_dir, "lipid_annotated_feature_matrix_ms1_ms2_comparespectra.csv"),
  row.names = FALSE
)

write.csv(
  combined[combined$annotation_support_type == "MS1+MS2", , drop = FALSE],
  file = file.path(output_dir, "lipid_annotated_feature_matrix_ms1_ms2_comparespectra_supported_only.csv"),
  row.names = FALSE
)

write.csv(
  ms1_all,
  file = file.path(output_dir, "lipid_annotated_feature_matrix_ms1_candidates_long.csv"),
  row.names = FALSE
)

write.csv(
  ms1_feature_summary,
  file = file.path(output_dir, "lipid_annotated_feature_matrix_ms1_candidates_summary.csv"),
  row.names = FALSE
)

write.csv(
  summary_table,
  file = file.path(output_dir, "lipid_annotated_feature_matrix_ms1_ms2_comparespectra_summary.csv"),
  row.names = FALSE
)

print(paste("MS1 annotated features:", summary_table$ms1_annotated_features))
print(paste("MS2 annotated features:", summary_table$ms2_annotated_features))
print(paste("MS1+MS2 supported features:", summary_table$ms1_ms2_annotated_features))
```

## 9. PCA と loadings でサンプル差と寄与脂質を読む

統合した annotation matrix を使って PCA を行い、
サンプル差とその差に効いている脂質候補を読み解きます。

### 入力ファイル・出力先・基本設定

``` r
rm(list = ls(all = TRUE))

if (!requireNamespace("loadings", quietly = TRUE)) {
  stop("Package 'loadings' is required.")
}

input_file <- "./output/lipid_annotation_matrix_combined_comparespectra/lipid_annotated_feature_matrix_ms1_ms2_comparespectra.csv"
output_dir <- "./output/lipid_annotation_matrix_combined_comparespectra/pca_loadings"

if (!file.exists(input_file)) stop("Annotated matrix was not found.")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

extract_algae_type <- function(sample_names) {
  x <- sub("^Posi_Ida_", "", sample_names)
  x <- sub("_[0-9]+$", "", x)
  x
}

annotation_cols <- c(
  "feature_id", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks",
  "isotopes", "adduct", "pcgroup", "is_monoisotope_like", "is_major_lipid_adduct",
  "ms1_candidate_count", "ms1_unique_compound_count", "ms1_query_adducts",
  "ms1_candidate_compound_ids", "ms1_candidate_names", "ms1_candidate_formulas",
  "ms1_category_codes", "ms1_min_abs_mz_diff", "ms1_min_abs_ppm_error",
  "ms2_library_entry_id", "ms2_annotation_name", "ms2_annotation_precursor_type",
  "ms2_annotation_compound_class", "ms2_sample_support", "ms2_mean_score",
  "ms2_max_score", "ms2_mean_matched_peak_count", "ms2_max_matched_peak_count",
  "annotation_support_type", "final_annotation_name", "final_annotation_source",
  "mzml"
)

mat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
sample_cols <- grep("^Posi_Ida_", names(mat), value = TRUE)
if (length(sample_cols) < 2) stop("At least two sample columns are required for PCA.")

feature_annotations <- mat[, intersect(annotation_cols, names(mat)), drop = FALSE]
rownames(feature_annotations) <- mat$feature_id
feature_annotations <- feature_annotations[, setdiff(names(feature_annotations), "feature_id"), drop = FALSE]

intensity_matrix <- as.matrix(mat[, sample_cols, drop = FALSE])
rownames(intensity_matrix) <- mat$feature_id
mode(intensity_matrix) <- "numeric"
```

### 強度行列から分散がゼロの feature を除き、log 変換して PCA 用に整形する

``` r
# 1. 強度行列から分散がゼロの feature を除き、log 変換して PCA 用に整形する
nonzero_sd <- apply(intensity_matrix, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
intensity_matrix <- intensity_matrix[nonzero_sd, , drop = FALSE]
feature_annotations <- feature_annotations[rownames(intensity_matrix), , drop = FALSE]

log_matrix <- log1p(intensity_matrix)
pca_input <- t(log_matrix)
```

### PCA を計算してスコアと loading を取り出す

``` r
# 2. PCA を計算してスコアと loading を取り出す
pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
pca_result <- loadings::pca_loading(pca_result)

pca_scores <- as.data.frame(pca_result$x)
pca_scores$sample <- rownames(pca_scores)
pca_scores$algae_type <- extract_algae_type(pca_scores$sample)

loading_r <- as.data.frame(pca_result$loading$R)
loading_p <- as.data.frame(pca_result$loading$p.value)
names(loading_r) <- paste0(names(loading_r), "_loading")
names(loading_p) <- paste0(gsub("\\.loading$", "", names(loading_r)), "_p_value")

feature_loadings <- cbind(
  data.frame(feature_id = rownames(feature_annotations), stringsAsFactors = FALSE),
  feature_annotations,
  loading_r,
  loading_p
)

importance <- as.data.frame(summary(pca_result)$importance)
importance$metric <- rownames(importance)
rownames(importance) <- NULL

top_pc1 <- feature_loadings[order(-abs(feature_loadings$PC1_loading)), , drop = FALSE]
top_pc2 <- feature_loadings[order(-abs(feature_loadings$PC2_loading)), , drop = FALSE]
```

### 結果をファイルに保存する

``` r
# 3. 結果をファイルに保存する
write.csv(
  pca_scores,
  file = file.path(output_dir, "lipid_annotation_pca_scores.csv"),
  row.names = FALSE
)
write.csv(
  feature_loadings,
  file = file.path(output_dir, "lipid_annotation_pca_feature_loadings.csv"),
  row.names = FALSE
)
write.csv(
  head(top_pc1, 100),
  file = file.path(output_dir, "lipid_annotation_pca_top100_pc1_loadings.csv"),
  row.names = FALSE
)
write.csv(
  head(top_pc2, 100),
  file = file.path(output_dir, "lipid_annotation_pca_top100_pc2_loadings.csv"),
  row.names = FALSE
)
write.csv(
  importance,
  file = file.path(output_dir, "lipid_annotation_pca_importance.csv"),
  row.names = FALSE
)
saveRDS(
  pca_result,
  file = file.path(output_dir, "lipid_annotation_pca_result.rds")
)
```

### PCA スコアプロットを PNG で保存する

``` r
# 4. PCA スコアプロットを PNG で保存する
algae_types <- unique(pca_scores$algae_type)
palette_colors <- setNames(grDevices::hcl.colors(length(algae_types), "Dark 3"), algae_types)
pc1_var <- round(100 * summary(pca_result)$importance[2, 1], 1)
pc2_var <- round(100 * summary(pca_result)$importance[2, 2], 1)

grDevices::png(
  filename = file.path(output_dir, "lipid_annotation_pca_score_plot.png"),
  width = 1800,
  height = 1400,
  res = 200
)
graphics::par(mar = c(5, 5, 3, 2))
graphics::plot(
  pca_scores$PC1,
  pca_scores$PC2,
  col = palette_colors[pca_scores$algae_type],
  pch = 19,
  cex = 1.4,
  xlab = paste0("PC1 (", pc1_var, "%)"),
  ylab = paste0("PC2 (", pc2_var, "%)"),
  main = "PCA Score Plot"
)
graphics::text(
  pca_scores$PC1,
  pca_scores$PC2,
  labels = pca_scores$sample,
  pos = 3,
  cex = 0.45
)
grDevices::dev.off()

cat("Samples:", nrow(pca_input), "\n")
cat("Features used for PCA:", ncol(pca_input), "\n")
```
