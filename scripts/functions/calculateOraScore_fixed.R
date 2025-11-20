calculateOraScore_fixed <- function(mSetObj = NA,
                                     nodeImp = "rbc",
                                     method = "hyperg") {

  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)

  if (!exists("current.kegglib")) {
    stop("❌ KEGG library not loaded. Run DownloadKEGGoffline() first.")
  }

  # -------------------------------
  # 1. mapping 後的 metabolite list
  # -------------------------------
  nm.map <- GetFinalNameMap(mSetObj)

  if (mSetObj$pathwaylibtype == "KEGG") {
    valid.inx <- !(is.na(nm.map$kegg) | duplicated(nm.map$kegg))
    ora.vec <- nm.map$kegg[valid.inx]
  } else {
    stop("❌ Only KEGG supported in bypass mode.")
  }

  q.size <- length(ora.vec)
  if (q.size < 2) stop("❌ Not enough KEGG IDs.")

  # -------------------------------
  # 2. load offline KEGG library
  # -------------------------------
  current.mset <- current.kegglib$mset.list
  uniq.count <- current.kegglib$uniq.count

  # -------------------------------
  # 3. filter universe
  # -------------------------------
  my.univ <- unique(unlist(current.mset))
  ora.vec <- ora.vec[ora.vec %in% my.univ]
  q.size <- length(ora.vec)

  if (q.size < 3) stop("❌ Too few matched metabolites for ORA.")

  # -------------------------------
  # 4. calculate hits
  # -------------------------------
  hits <- lapply(current.mset, function(x) x[x %in% ora.vec])
  hit.num <- sapply(hits, length)
  set.num <- sapply(current.mset, length)

  # -------------------------------
  # 5. build result matrix
  # -------------------------------
  res.mat <- matrix(0, nrow = length(current.mset), ncol = 8)
  colnames(res.mat) <- c("Total","Expected","Hits","Raw p",
                         "-log10(p)","Holm","FDR","Impact")
  rownames(res.mat) <- names(current.mset)

  # Topology score
  imp.list <- if (nodeImp == "rbc") current.kegglib$rbc else current.kegglib$dgr
  imp.list <- imp.list[names(current.mset)]

  res.mat[,1] <- set.num
  res.mat[,2] <- q.size * (set.num / uniq.count)
  res.mat[,3] <- hit.num

  # Raw P
  if (method == "fisher") {
    res.mat[,4] <- mapply(function(k, K, n, N) {
      fisher.test(matrix(c(k, K-k, n-k, N-K-(n-k)), 2))$p.value
    }, hit.num, set.num, q.size, uniq.count)
  } else {
    res.mat[,4] <- phyper(hit.num - 1, set.num,
                          uniq.count - set.num,
                          q.size, lower.tail = FALSE)
  }

  res.mat[,5] <- -log10(res.mat[,4])
  res.mat[,6] <- p.adjust(res.mat[,4], method = "holm")
  res.mat[,7] <- p.adjust(res.mat[,4], method = "fdr")
  res.mat[,8] <- mapply(function(x,y) sum(x[y]), imp.list, hits)

  res.mat <- res.mat[hit.num > 0, ]
  ord <- order(res.mat[,4], res.mat[,8])
  res.mat <- res.mat[ord, ]

  # save back
  mSetObj$analSet$ora.mat <- signif(res.mat, 5)
  mSetObj$analSet$ora.hits <- hits
  mSetObj$analSet$node.imp <- nodeImp

  return(MetaboAnalystR:::.set.mSet(mSetObj))
}
