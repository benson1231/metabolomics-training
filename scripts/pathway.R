# This example data set is a peak list from a COVID-19 study
# Two groups (COVID-19 vs. Healthy Controls) are included for this example

download.file("https://raw.githubusercontent.com/Zhiqiang-PANG/MetaboRaw/master/examples/peaks_ms1.txt", 
              destfile = "peaks_ms1.txt", 
              mode = "auto")

# Load MetaboAnalystR
suppressPackageStartupMessages({
    library(MetaboAnalystR)
})

# Create objects for storing processed data
mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
# Set Peak format as "mpt" here.
# This option can be "mp", "mt", "mpt", "mprt", "mrt", "mpr", "rmp", "rmt".
# "rmp" and "rmt" refers to peaks ranked by p values or t scores;
# For other options, "m" means "mz"; "p" means "p value"; "t" means "t score"; "r" means "retention time";
mSet <- SetPeakFormat(mSet, "mpt")

# Set parameters for analysis, in this case the mass accuracy is set to 15 ppm, 
# the mode of the MS instrument is "mixed" (contains both positive and negative), 
# and the p-value cut-off is 0.02
mSet <- UpdateInstrumentParameters(mSet, 15.0, "mixed", "yes", 0.02)

# Read in peak-list data
mSet <- Read.PeakListData(mSet, "peaks_ms1.txt");

# set retention time included, unit is "seconds"
mSet <- SetRTincluded(mSet, "seconds")

# Sanity check of the uploaded data
mSet <- SanityCheckMummichogData(mSet)

# set peak enrichment method. Method can be one of the "mum", "gsea" or "integ";
# Method, "integ" means integration of both mummichog and GSEA algorithm;
# version can be "v1" or "v2" ("v1" will use m/z only; "v2" will use both "m/z" and "retention time")
mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")

# Here we use the top 10% peaks as the p value cutoff
pval <- sort(mSet[["dataSet"]][["mummi.proc"]][["p.value"]])[ceiling(length(mSet[["dataSet"]][["mummi.proc"]][["p.value"]])*0.1)]
mSet <- SetMummichogPval(mSet, pval)

# Perform the mummichog algorith, in this case the model is the human MFN model, using "current" version by default
# This function may take sometime for processing, and will output the pathway-results and the compound matching tables in your working directory
mSet <- PerformPSEA(mSet, "hsa_mfn", "current", 3 , 100)

# To view the results of the pathway analysis
mSet <- PlotPeaks2Paths(mSet, "peaks_to_paths_ms1_", "png", 72, width=8)