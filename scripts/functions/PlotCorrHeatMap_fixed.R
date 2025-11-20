PlotCorrHeatMap_fixed <- function(
  mSetObj,
  out_png,
  target      = "col",        # "row" = feature correlation, "col" = sample correlation
  method      = "pearson",    # pearson / spearman
  corrCutoff  = 0,            # threshold for filtering
  width       = 10,
  height      = 10,
  fontsize    = 10
){
  if(!requireNamespace("pheatmap", quietly = TRUE)){
    stop("Please install pheatmap: install.packages('pheatmap')")
  }
  
  # -----------------------------
  # 1. 取 normalized matrix
  # -----------------------------
  if(!"norm" %in% names(mSetObj$dataSet)){
    stop("mSetObj$dataSet$norm is missing — ensure you ran normalization.")
  }
  
  mat <- mSetObj$dataSet$norm
  
  # -----------------------------
  # 2. row/col 模式（feature/samples）
  # -----------------------------
  if(target == "row"){
    mat <- t(mat)
  }
  
  # -----------------------------
  # 3. 計算 correlation matrix
  # -----------------------------
  corr.mat <- cor(mat, method = method)
  
  # -----------------------------
  # 4. correlation cutoff
  # -----------------------------
  corr.mat[abs(corr.mat) < corrCutoff] <- 0
  
  # -----------------------------
  # 5. 顏色（與 MetaboAnalyst 相似）
  # -----------------------------
  colors <- colorRampPalette(
    c("#0571b0","#92c5de","white","#f4a582","#ca0020")
  )(256)
  
  # -----------------------------
  # 6. 產生 PNG heatmap
  # -----------------------------
  pheatmap::pheatmap(
    corr.mat,
    filename      = out_png,
    color         = colors,
    cluster_rows  = TRUE,
    cluster_cols  = TRUE,
    fontsize      = fontsize,
    width         = width,
    height        = height
  )
  
  message("✔ Correlation heatmap saved to: ", out_png)
  return(mSetObj)
}
