# ---------------------------------------------------------------
# Pathway Analysis using Mummichog (export → pathway/)
# ---------------------------------------------------------------

data_dir <- "results/pathway"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Start log (terminal + file)
sink("results/pathway/pathway.log", split = TRUE)

# Download example MS1 peak list
download.file(
  "https://raw.githubusercontent.com/Zhiqiang-PANG/MetaboRaw/master/examples/peaks_ms1.txt",
  destfile = file.path(data_dir, "peaks_ms1.txt"),
  mode = "auto"
)

suppressPackageStartupMessages(library(MetaboAnalystR))

# Helper function to redirect all outputs to pathway/
redirect_path <- function(filename) file.path(data_dir, filename)

# ---------------------------------------------------------------
# Initialize MetaboAnalyst data object
# ---------------------------------------------------------------

mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
mSet <- SetPeakFormat(mSet, "mpt")
mSet <- UpdateInstrumentParameters(mSet, 15.0, "mixed", "yes", 0.02)

# Read peak list
mSet <- Read.PeakListData(mSet, "peaks_ms1.txt")

# Set RT unit
mSet <- SetRTincluded(mSet, "seconds")

# QC / sanity check
mSet <- SanityCheckMummichogData(mSet)

# Mummichog enrichment settings
mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")

# Top 10% p-value cutoff
pval_vec <- mSet$dataSet$mummi.proc$p.value
pval_cut <- sort(pval_vec)[ceiling(length(pval_vec) * 0.1)]
mSet <- SetMummichogPval(mSet, pval_cut)

# ---------------------------------------------------------------
# Perform Mummichog pathway enrichment
# ---------------------------------------------------------------
mSet <- PerformPSEA(mSet, "hsa_mfn", "current", 3 , 100)


# ---------------------------------------------------------------
# Export pathway visualization (saved to pathway/)
# ---------------------------------------------------------------

mSet <- PlotPeaks2Paths(mSet, "peaks_to_paths_ms1_", "png", 72, width=8)

message("✔ All pathway analysis outputs saved to: ", data_dir)

# Stop logging
sink()
