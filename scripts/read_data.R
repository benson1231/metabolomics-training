### Tutorial from https://www.metaboanalyst.ca/resources/vignettes/Introductions.html
suppressPackageStartupMessages(library(MetaboAnalystR))

# Create output folder
dir.create("results/read_data", recursive = TRUE, showWarnings = FALSE)

# Start log (capture ALL outputs: stdout + messages)
log_file <- "results/read_data/read_data.log"
log_con <- file(log_file, open = "wt")

sink(log_con, type = "output")    # capture printed output
sink(log_con, type = "message")   # capture warnings/messages

cat("====== MetaboAnalystR Read Data Log ======\n")
cat("Start time:", as.character(Sys.time()), "\n\n")

###################################################
### 1. Read compound concentration dataset (conc)
###################################################

cat(">>> Step 1: Reading compound concentration dataset...\n")

mSet <- InitDataObjects("conc", "stat", FALSE)

conc_file <- "results/read_data/human_cachexia.csv"
download.file(
  "https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv",
  destfile = conc_file
)

mSet <- Read.TextData(mSet, conc_file, "rowu", "disc")

cat("\n--- Read Message ---\n")
print(mSet$msgSet$read.msg)
cat("\n")

cat("✓ Finished reading concentration dataset.\n\n")

##########################################
### 2. Read peak intensity table
##########################################

cat(">>> Step 2: Reading LC-MS peak intensity table...\n")

mSet <- InitDataObjects("pktable", "stat", FALSE)

peak_file <- "results/read_data/lcms_table.csv"
download.file(
  "https://rest.xialab.ca/api/download/metaboanalyst/lcms_table.csv",
  destfile = peak_file
)

mSet <- Read.TextData(mSet, peak_file, "rowu", "disc")

cat("\n--- Read Message ---\n")
print(mSet$msgSet$read.msg)
cat("\n")

cat("✓ Finished reading LC-MS table.\n\n")

##########################################
### 3. Import example 3-column MS peaks
##########################################

cat(">>> Step 3: Reading example LC-MS 3-column peaks...\n")

# Clean upload folder (do NOT clean workspace)
unlink("results/read_data/upload", recursive = TRUE)

mSet <- InitDataObjects("mspeak", "stat", FALSE)

zip_path <- "results/read_data/lcms_3col_peaks.zip"
download.file(
  "https://rest.xialab.ca/api/download/metaboanalyst/lcms_3col_peaks.zip",
  destfile = zip_path
)

UnzipUploadedFile(zip_path, "results/read_data/upload", FALSE)
mSet <- Read.PeakList(mSet, "results/read_data/upload")

cat("✓ Finished importing MS peak list.\n\n")

###################################################
### End of log
###################################################

cat("End time:", as.character(Sys.time()), "\n")
cat("====== Log End ======\n")

sink(type = "message")
sink(type = "output")
close(log_con)
