### https://www.metaboanalyst.ca/resources/vignettes/Pathway_Analysis.html

### Pathway Analysis — Clean Working Example
rm(list = ls())

# -----------------------------------------
# Setup directories
# -----------------------------------------
data_dir <- "results/pathway"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(data_dir, "*"), recursive = TRUE)

# -----------------------------------------
# Load MetaboAnalystR
# -----------------------------------------
suppressPackageStartupMessages(library(MetaboAnalystR))

# -----------------------------------------
# Compound list for ORA / Pathway Analysis
# -----------------------------------------
tmp.vec <- c(
  "Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine",
  "Fumaric acid", "Glycine", "Homocysteine", "L-Cysteine",
  "L-Isolucine", "L-Phenylalanine", "L-Serine", "L-Threonine",
  "L-Tyrosine", "L-Valine", "Phenylpyruvic acid", "Propionic acid",
  "Pyruvic acid", "Sarcosine"
)

# -----------------------------------------
# Initialize for pathway ORA
# -----------------------------------------
mSet <- InitDataObjects("conc", "pathora", FALSE)

# Load compound list
mSet <- Setup.MapData(mSet, tmp.vec)

# Cross-reference names vs HMDB/PubChem/KEGG/etc.
mSet <- CrossReferencing(mSet, "name")

# Create mapping results table
mSet <- CreateMappingResultTable(mSet)

# -----------------------------------------
# Choose KEGG pathway library
# hsa = human
# -----------------------------------------
mSet <- SetKEGG.PathLib(mSet, "hsa", "current")

# Disable metabolome filtering (include all)
mSet <- SetMetabolomeFilter(mSet, FALSE)









# still in development -------------------------------------------------------------------
source("scripts/functions/DownloadKEGGoffline.R")
DownloadKEGGoffline()

source("scripts/functions/calculateOraScore_fixed.R")
mSet <- calculateOraScore_fixed(mSet, "rbc", "hyperg")

### API error -----------------------------------------------------------------------------
# -----------------------------------------
# Calculate ORA Hypergeometric Test
# -----------------------------------------
### API error
mSet <- CalculateOraScore(mSet, "rbc", "hyperg")

# -----------------------------------------
# Plot Pathway Overview
# -----------------------------------------
mSet<-PlotPathSummary(mSet, T, file.path(data_dir, "path_view_0_"), "png", 72, width=NA)

message("✔ Pathway overview saved.")

# -----------------------------------------
# Plot KEGG Pathway Diagram (選擇特定路徑)
# Example: Glycine, serine and threonine metabolism
# -----------------------------------------
mSet<-PlotKEGGPath(mSet, "Glycine, serine and threonine metabolism",576, 480, "png", NULL)

message("✔ Specific KEGG pathway saved.")
