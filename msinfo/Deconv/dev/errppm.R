# --------------------
#   エラーppm 関数
# --------------------
e_ppm <- function(mz_m, mz_t){
  # mz_m:測定, mz_t:理論
  # e_ppm <- abs ((mz_m-mz_t)/mz_m)*10^6
  e_ppm <- ((mz_m-mz_t)/mz_t)*10^6
  return(e_ppm)
}

#y <- e_ppm(x1,mz) # 目的変数




44.998201