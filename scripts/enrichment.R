### https://www.metaboanalyst.ca/resources/vignettes/Enrichment_Analysis.html

rm(list = ls())

suppressPackageStartupMessages(library(MetaboAnalystR))

data_dir <- "results/enrichment"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(data_dir, "*"), recursive = TRUE)

## When input is a list
# Create vector consisting of compounds for enrichment analysis 
tmp.vec <- c("Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine", "Fumaric acid", "Glycine", "Homocysteine", "L-Cysteine", "L-Isolucine", "L-Phenylalanine", "L-Serine", "L-Threonine", "L-Tyrosine", "L-Valine", "Phenylpyruvic acid", "Propionic acid", "Pyruvic acid", "Sarcosine", "Arsenic", "Benzene", "Caffeic acid", "Cotinine", "Cadmium", "Lead", "Thiocyanate")

mSet <- InitDataObjects("conc", "msetora", FALSE)
mSet <- Setup.MapData(mSet, tmp.vec)
mSet <- CrossReferencing(mSet, "name")
mSet <- CreateMappingResultTable(mSet)
mSet <- SetMetabolomeFilter(mSet, F)
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 0)
mSet <- CalculateHyperScore(mSet)
mSet <- PlotORA(mSet, file.path(data_dir, "ora_0_"), "bar", "png", 72)





### Quantitative Enrichment Analysis
suppressPackageStartupMessages(library(MetaboAnalystR))

rm(list = ls())

data_dir <- "results/enrichment"

#------------------------------------------
#  1. 初始化 QEA
#------------------------------------------
mSet <- InitDataObjects("conc", "msetqea", FALSE)

# 載入 cachexia 濃度表（MetaboAnalyst 官方示範資料）
conc_file <- file.path(data_dir, "human_cachexia.csv")
download.file(
  "https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv",
  destfile = conc_file,
  mode = "wb"
)

mSet <- Read.TextData(mSet, conc_file, "rowu", "disc")

#------------------------------------------
#  2. 交叉比對 metabolite name
#------------------------------------------
mSet <- CrossReferencing(mSet, "name")
mSet <- CreateMappingResultTable(mSet)

#------------------------------------------
#  3. 資料檢查、補 missing
#------------------------------------------
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)

#------------------------------------------
#  4. 不做 normalization（QEA 免）
#------------------------------------------
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "NULL", "NULL", NULL, ratio=FALSE)

#------------------------------------------
#  5. 選 pathway library
#------------------------------------------
mSet <- SetMetabolomeFilter(mSet, FALSE)
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 0)

#------------------------------------------
#  6. 計算 QEA
#------------------------------------------
mSet <- CalculateGlobalTestScore(mSet)

#------------------------------------------
#  7. 畫 QEA overview
#------------------------------------------
mSet <- PlotQEA.Overview(
  mSet,
  imgName = file.path(data_dir, "qea_0_"),
  format = "png",
  dpi = 72,
  width = NA
)

message("✔ QEA finished and saved to: ", data_dir)