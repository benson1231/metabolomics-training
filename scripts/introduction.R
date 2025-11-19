### Tutorial from https://www.metaboanalyst.ca/resources/vignettes/Introductions.html
# Load MetaboAnalystR
library(MetaboAnalystR)

# Create output folder
dir.create("results/introduction", recursive = TRUE, showWarnings = FALSE)

# Start log (terminal + file)
sink("results/introduction/metabo_introduction.log", split = TRUE)

cat("====== MetaboAnalystR Introduction Log ======\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

###################################################
### 1. Initialize and download NMR peak dataset
###################################################
cat(">>> Step 1: Downloading NMR peak dataset...\n")

mSet <- InitDataObjects("nmrpeak", "stat", FALSE)

zip_file <- "results/introduction/nmr_peaks.zip"
download.file(
    "https://rest.xialab.ca/api/download/metaboanalyst/nmr_peaks.zip",
    destfile = zip_file,
    method = "auto"
)

cat("✓ Download completed.\n\n")

###################################################
### 2. Unzip input files
###################################################
cat(">>> Step 2: Unzipping dataset...\n")

upload_dir <- "results/introduction/upload"
UnzipUploadedFile(zip_file, upload_dir)

cat("✓ Unzip completed.\n\n")

###################################################
### 3. Read data into mSet
###################################################
cat(">>> Step 3: Reading peak list...\n")

mSet <- Read.PeakList(mSet, upload_dir)
cat("✓ Read completed.\n\n")


###################################################
### 4. Peak grouping
###################################################
cat(">>> Step 4: Peak grouping...\n")

mSet <- GroupPeakList(mSet, 0.025, 30.0)
mSet <- SetPeakList.GroupValues(mSet)

cat("✓ Peak grouping completed.\n\n")


###################################################
### 5. Sanity check
###################################################
cat(">>> Step 5: Running sanity check...\n")

mSet <- SanityCheckData(mSet)
cat("✓ Sanity check completed.\n\n")


###################################################
### 6. Missing value handling
###################################################
cat(">>> Step 6: Handling missing values...\n")

mSet <- ReplaceMin(mSet)

writeLines(
    mSet$msgSet$replace.msg,
    "results/introduction/replace_min_messages.txt"
)

cat("✓ Missing value replacement completed.\n\n")


###################################################
### 7. KNN imputation
###################################################
cat(">>> Step 7: Performing KNN imputation...\n")

mSet <- ImputeMissingVar(mSet, method = "knn_smp")

cat("✓ KNN imputation completed.\n\n")


###################################################
### 8. Small sample size check
###################################################
cat(">>> Step 8: Checking sample size...\n")

small_check <- IsSmallSmplSize(mSet)
writeLines(
    paste("Small sample size check:", small_check),
    "results/introduction/small_sample_check.txt"
)

cat("✓ Sample size check completed.\n\n")


###################################################
### 9. Normalization
###################################################
cat(">>> Step 9: Performing normalization (Quantile + Log + Mean Center)...\n")

mSet <- PreparePrenormData(mSet)

mSet<-Normalization(mSet, "SamplePQN", "NULL", "NULL", "P037", ratio=FALSE, ratioNum=20)

cat("✓ Normalization completed.\n\n")


###################################################
### 10. Output normalized plots
###################################################
cat(">>> Step 10: Exporting normalization plots...\n")

mSet <- PlotNormSummary(
    mSet,
    "results/introduction/feature_norm.png",
    format = "png",
    dpi = 72
)

mSet <- PlotSampleNormSummary(
    mSet,
    "results/introduction/sample_norm.pdf",
    format = "pdf"
)

cat("✓ Plots saved.\n\n")


###################################################
### 11. Variable filtering examples
###################################################
cat(">>> Step 11: Filtering variables...\n")

mSet <- FilterVariable(mSet, "mad", 5, "F", 25, TRUE)
mSet <- FilterVariable(mSet, "nrsd", 5, "T", 25, TRUE)

cat("✓ Filtering completed.\n\n")


###################################################
### 12. Update dataset
###################################################
cat(">>> Step 12: Updating dataset...\n")

# Remove a sample from the data set, in this case sample "PIF_178"
smpl.nm.vec <- c("PIF_178")# used to remove certain samples

# Remove a feature from the data set
feature.nm.vec <- c("2-Aminobutyrate")# used to remove certain feature, i.e. 2-Aminobutyrate

# Remove a group from the data set, in this case remove the "control" samples
grp.nm.vec <- c("control") # used to retain certain groups

mSet <- UpdateData(mSet)

cat("✓ Dataset updated.\n\n")


###################################################
### 13. Save transformed data
###################################################
cat(">>> Step 13: Saving transformed data...\n")

setwd("results/introduction")
PreparePDFReport(mSet, "User Name")

SaveTransformedData(mSet)

cat("✓ Transformed data saved.\n\n")

###################################################
### End of pipeline
###################################################
cat("End time:", as.character(Sys.time()), "\n")
cat("====== Log End ======\n")

# Stop logging
sink()
