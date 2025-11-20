DownloadKEGGoffline <- function() {

  message("📦 Downloading KEGG offline libraries...")

  base_url <- "https://raw.githubusercontent.com/xialab/MetaboAnalystR/master/inst/extdata/libs/"

  libs <- c(
    "kegg_hsa.qs",
    "kegg_mmu.qs",
    "kegg_pathway.qs",
    "kegg_graph.qs",
    "kegg_network.qs"
  )

  dest_dir <- file.path(find.package("MetaboAnalystR"), "extdata", "libs")
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  for (f in libs) {
    url <- paste0(base_url, f)
    dest <- file.path(dest_dir, f)

    message("⬇  Downloading: ", f)

    tryCatch({
      download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
    }, error = function(e) {
      message("❌ Failed: ", f)
    })
  }

  message("✔ KEGG offline libraries installed at: ", dest_dir)
}
