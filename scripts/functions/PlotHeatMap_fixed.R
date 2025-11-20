# --------------------------------------------------
# identical to MetaboAnalystR scale_mat()
# --------------------------------------------------
scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
} 
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

PlotHeatMap_fixed <- function(
    mSetObj,
    imgName,
    dataOpt,
    scaleOpt,
    smplDist,
    clstDist,
    palette,
    grp.ave       = FALSE,
    var.inx       = NULL,
    rowV          = TRUE,
    colV          = TRUE,
    showColnm     = TRUE,
    showRownm     = TRUE,
    width         = 10,
    height        = 12,
    fontsize      = 8
){

    if(!requireNamespace("pheatmap", quietly = TRUE)){
        stop("Please install pheatmap: install.packages('pheatmap')")
    }

    # -----------------------------
    # 1. Load dataset
    # -----------------------------
    if(dataOpt == "norm"){
        dat <- mSetObj$dataSet$norm
    } else {
        dat <- qs::qread("prenorm.qs")
    }

    # 選 feature
    if(!is.null(var.inx)){
        dat <- dat[, var.inx, drop = FALSE]
    }

    cls <- mSetObj$dataSet$cls

    # -----------------------------
    # 2. Tag annotation（與原生流程一致）
    # -----------------------------
    annotation <- data.frame(class = cls)
    rownames(annotation) <- rownames(dat)

    # -----------------------------
    # 3. Group average (原邏輯完全一致)
    # -----------------------------
    if(grp.ave){
        lv <- levels(cls)
        mat2 <- matrix(nrow = length(lv), ncol = ncol(dat))
        for(i in seq_along(lv)){
            mat2[i, ] <- colMeans(dat[cls == lv[i], , drop = FALSE])
        }
        dat <- mat2
        rownames(dat) <- lv
        annotation <- data.frame(class = factor(lv, levels = lv))
        rownames(annotation) <- lv
    }

    # -----------------------------
    # 4. Scale 與 transpose（與 PlotHeatMap 完全一致）
    # -----------------------------
    dat <- t(dat)
    dat <- scale_mat(dat, scaleOpt)
    dat <- round(dat, 5)

    # -----------------------------
    # 5. Heatmap 色系（與 MetaboAnalyst 一致）
    # -----------------------------
    if(palette == "bwm"){
        hm.colors <- colorRampPalette(
            c("#0571b0","#92c5de","white","#f4a582","#ca0020")
        )(256)
    } else {
        stop("目前先支援 palette='bwm'，需要其他我可幫你加")
    }

    # -----------------------------
    # 6. Distance / Clustering（與原流程一致）
    # -----------------------------
    if(smplDist == "correlation"){
        dist_row <- as.dist(1 - cor(dat, method="pearson"))
        dist_col <- as.dist(1 - cor(t(dat), method="pearson"))
    } else {
        dist_row <- dist(dat, method = smplDist)
        dist_col <- dist(t(dat), method = smplDist)
    }

    clust_row <- if(rowV) hclust(dist_row, method = clstDist) else FALSE
    clust_col <- if(colV) hclust(dist_col, method = clstDist) else FALSE

    # -----------------------------
    # 7. 固定色階 -4 ~ 4
    # -----------------------------
    breaks <- seq(-4, 4, length.out = 257)

    # -----------------------------
    # 8. Export PNG
    # -----------------------------
    pheatmap::pheatmap(
        dat,
        filename      = imgName,
        color         = hm.colors,
        breaks        = breaks,
        annotation_col = annotation,
        cluster_rows  = clust_row,
        cluster_cols  = clust_col,
        show_colnames = showColnm,
        show_rownames = showRownm,
        fontsize      = fontsize,
        width         = width,
        height        = height,
        border_color  = "grey80"
    )

    message("✔ Static heatmap saved to: ", imgName)
    return(mSetObj)
}
