### Tutorial from https://www.metaboanalyst.ca/resources/vignettes/LCMSMS_Raw_Spectral_Processing.html
### LC-MS/MS Raw Spectral Processing Pipeline
### Fully annotated + Logging enabled for reproducibility

suppressPackageStartupMessages({
    library(MetaboAnalystR)
    library(OptiLCMS)
    library(ggplot2)
})

# ---------------------------------------------------------------
# 0. Global Setup
# ---------------------------------------------------------------

data_dir <- "results/lcms"
options(timeout = 600)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------
# Enable Logging (sink)
# ---------------------------------------------------------------

export_dir <- file.path(data_dir, "exports")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(data_dir, "lcms.log")
log_con <- file(log_file, open = "wt")

sink(log_con, type = "output")
sink(log_con, type = "message")

log_timestamp <- function(msg) {
  cat("\n===============================\n")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n")
  cat("===============================\n\n")
}

log_timestamp("Pipeline started")


# ---------------------------------------------------------------
# 1. Download example dataset
# ---------------------------------------------------------------
log_timestamp("Step 1: Downloading malaria_raw.zip")

download.file(
  "https://api2.xialab.ca/api/download/metaboanalyst/malaria_raw.zip",
  destfile = file.path(data_dir, "malaria_raw.zip"),
  mode = "wb"
)


# ---------------------------------------------------------------
# 2. Extract into upload/
# ---------------------------------------------------------------
log_timestamp("Step 2: Unzipping raw data")

upload_dir <- file.path(data_dir, "upload")
unzip(file.path(data_dir, "malaria_raw.zip"), exdir = upload_dir)

mal_dir <- file.path(upload_dir, "malaria_LCMS")


# ---------------------------------------------------------------
# 3. Create group folders
# ---------------------------------------------------------------
log_timestamp("Step 3: Creating group folders")

grp_naive_dir <- file.path(upload_dir, "Naive")
grp_qc_dir    <- file.path(upload_dir, "QC")
grp_semi_dir  <- file.path(upload_dir, "Semi_immune")

dir.create(grp_naive_dir, showWarnings = FALSE)
dir.create(grp_qc_dir,    showWarnings = FALSE)
dir.create(grp_semi_dir,  showWarnings = FALSE)


# ---------------------------------------------------------------
# 4. Distribute raw .zip into groups
# ---------------------------------------------------------------
log_timestamp("Step 4: Extracting mzML into group folders")

all_zips <- list.files(mal_dir, pattern = "\\.zip$", full.names = TRUE)

qc_zips    <- grep("^QC_",    basename(all_zips), value = TRUE)
naive_zips <- grep("^Naive_", basename(all_zips), value = TRUE)
semi_zips  <- grep("^Semi_",  basename(all_zips), value = TRUE)

for (zf in file.path(mal_dir, qc_zips)) unzip(zf, exdir = grp_qc_dir)
for (zf in file.path(mal_dir, naive_zips)) unzip(zf, exdir = grp_naive_dir)
for (zf in file.path(mal_dir, semi_zips))  unzip(zf, exdir = grp_semi_dir)

log_timestamp("All mzML files extracted")


# ---------------------------------------------------------------
# 5. ROI extraction (QC only)
# ---------------------------------------------------------------
log_timestamp("Step 5: ROI extraction using QC files")

qc_files <- list.files(grp_qc_dir, full.names = TRUE)
print(qc_files)

mSet <- PerformROIExtraction(qc_files, rt.idx = 0.9, rmConts = TRUE)


# ---------------------------------------------------------------
# 6. Parameter optimization
# ---------------------------------------------------------------
log_timestamp("Step 6: Parameter optimization")

best_params <- PerformParamsOptimization(
  mSet,
  param = SetPeakParam(platform = "UPLC-Q/E"),
  ncore = 4
)


# ---------------------------------------------------------------
# 7. Import all raw data
# ---------------------------------------------------------------
log_timestamp("Step 7: Importing ALL raw data")

mSet <- ImportRawMSData(
  path = upload_dir,
  plotSettings = SetPlotParam(Plot = TRUE)
)


# ---------------------------------------------------------------
# 8. Peak profiling (peak picking + alignment + filling)
# ---------------------------------------------------------------
log_timestamp("Step 8: Peak profiling")

mSet <- PerformPeakProfiling(
  mSet,
  Params = best_params,
  plotSettings = SetPlotParam(Plot = TRUE)
)


# ---------------------------------------------------------------
# 9. CAMERA annotation
# ---------------------------------------------------------------
log_timestamp("Step 9: CAMERA peak annotation")

annParams <- SetAnnotationParam(
  polarity = "positive",
  mz_abs_add = 0.015
)

mSet <- PerformPeakAnnotation(mSet, annParams)


# ---------------------------------------------------------------
# 10. Export results
# ---------------------------------------------------------------
log_timestamp("Step 10: Exporting results")

feat <- mSet@peakAnnotation$camera_output
write.csv(feat, file.path(export_dir, "camera_output_full.csv"), row.names = FALSE)

int_cols <- grep("\\.mzML$", colnames(feat))
metabo_input <- feat[, c("mz", "rt", int_cols)]
write.csv(metabo_input, file.path(export_dir, "metaboanalyst_input.csv"), row.names = FALSE)

if (!is.null(mSet@peakpicking$chromPeaks)) {
  write.csv(mSet@peakpicking$chromPeaks,
            file.path(export_dir, "chromPeaks.csv"),
            row.names = FALSE)
}

write.csv(mSet@peakAnnotation$groups,
          file.path(export_dir, "peak_groups.csv"),
          row.names = FALSE)

sink(file.path(export_dir, "pspectra_list.txt"))
print(mSet@peakAnnotation$AnnotateObject$pspectra)
sink()

write.csv(mSet@peakAnnotation$AnnotateObject$ruleset,
          file.path(export_dir, "annotation_ruleset.csv"),
          row.names = FALSE)

iso_list <- mSet@peakAnnotation$isoList
if (!is.null(iso_list)) {
  sink(file.path(export_dir, "isotope_list.txt"))
  print(iso_list)
  sink()
}

write.table(
  colnames(metabo_input)[3:ncol(metabo_input)],
  file.path(export_dir, "sample_names.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

rt_mz_summary <- feat[, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax")]
write.csv(rt_mz_summary, file.path(export_dir, "rt_mz_summary.csv"), row.names = FALSE)

write.csv(feat[, int_cols], file.path(export_dir, "intensity_matrix.csv"), row.names = FALSE)

saveRDS(mSet, file.path(export_dir, "mSet_processed.rds"))

# ---------------------------------------------------------------
# PCA Analysis (using CAMERA feature matrix)
# ---------------------------------------------------------------
log_timestamp("Step 11: PCA plotting")

# Extract intensity matrix (columns 7–18 in this dataset)
feat <- mSet@peakAnnotation$camera_output
int_mat <- feat[, 7:18]
int_mat[is.na(int_mat)] <- 0
int_log <- log2(int_mat + 1)

# Build metadata
sample_names <- colnames(int_log)
Group <- ifelse(grepl("Naive", sample_names, TRUE), "Naive",
         ifelse(grepl("QC", sample_names, TRUE), "QC", "Semi_immune"))

meta <- data.frame(Sample = sample_names, Group = factor(Group))

# PCA
pca_res <- prcomp(t(int_log), scale. = TRUE)
pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = meta$Group
)

# Plot
library(ggplot2)
p <- ggplot(pca_df, aes(PC1, PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.6, size = 3) +
  theme_bw() +
  ggtitle("PCA of LC-MS Feature Matrix (CAMERA output)")

ggsave(file.path(export_dir, "PCA.png"), plot = p, width = 8, height = 6, dpi = 300)

log_timestamp("PCA saved.")

# ---------------------------------------------------------------
# Finalize log
# ---------------------------------------------------------------
log_timestamp("Pipeline Completed Successfully")

sink(type = "output")
sink(type = "message")
close(log_con)
