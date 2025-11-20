PlotPCA3DScoreImg_fixed <- function(
  mSetObj = NA,
  imgName,
  format = "png",
  dpi = 72,
  width = NA,
  pc1 = 1,
  pc2 = 2,
  pc3 = 3,
  angle = 40
){

  # ← 使用 MetaboAnalystR internal function
  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)

  # Labels
  xlabel = paste0("PC", pc1, " (", round(100 * mSetObj$analSet$pca$variance[pc1], 1), "%)")
  ylabel = paste0("PC", pc2, " (", round(100 * mSetObj$analSet$pca$variance[pc2], 1), "%)")
  zlabel = paste0("PC", pc3, " (", round(100 * mSetObj$analSet$pca$variance[pc3], 1), "%)")

  imgName = paste0(imgName, "dpi", dpi, ".", format)

  # Width / height
  if (is.na(width)) {
    w <- 9
  } else if (width == 0) {
    w <- 7.2
  } else {
    w <- width
  }
  h <- w

  mSetObj$imgSet$pca.score3d <- imgName

  # ========== FIX COLOR BUG ==========
  cls <- mSetObj$dataSet$cls
  cls.factor <- as.factor(cls)
  palette <- grDevices::rainbow(length(levels(cls.factor)))

  if (!is.null(mSetObj$dataSet$cls.col)) {
    cols <- mSetObj$dataSet$cls.col
  } else if (!is.null(mSetObj$dataSet$cls.color)) {
    cols <- mSetObj$dataSet$cls.color
  } else {
    cols <- palette[cls.factor]
  }
  # ====================================

  pchs <- as.numeric(cls.factor)
  uniq.pchs <- unique(pchs)

  # Output device
  Cairo::Cairo(
    file = imgName,
    unit = "in",
    dpi = dpi,
    width = w,
    height = h,
    type = format,
    bg = "white"
  )

  # Plot
  MetaboAnalystR:::Plot3D(
    mSetObj$analSet$pca$x[, pc1],
    mSetObj$analSet$pca$x[, pc2],
    mSetObj$analSet$pca$x[, pc3],
    xlab = xlabel,
    ylab = ylabel,
    zlab = zlabel,
    angle = angle,
    color = cols,
    pch = pchs
  )

  legend("topleft", legend = levels(cls.factor), pch = uniq.pchs, col = palette)

  dev.off()

  return(MetaboAnalystR:::.set.mSet(mSetObj))
}
