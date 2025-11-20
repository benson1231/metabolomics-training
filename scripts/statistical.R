# Load MetaboAnalystR
suppressPackageStartupMessages(library(MetaboAnalystR))

# Clean global environment
rm(list = ls())

data_dir <- "results/statistical"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

library(MetaboAnalystR)
mSet<-InitDataObjects("conc", "stat", FALSE);

conc_file <- file.path(data_dir, "human_cachexia.csv")
download.file(
  "https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv",
  destfile = conc_file
)

mSet <- Read.TextData(mSet, conc_file, "rowu", "disc")
mSet <- SanityCheckData(mSet)


mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet);
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20);

mSet<-PlotNormSummary(mSet, file.path(data_dir, "norm_0_"), format ="png", dpi=72, width=NA);
mSet<-PlotSampleNormSummary(mSet, file.path(data_dir, "snorm_0_"), format = "png", dpi=72, width=NA);


### Fold-change analysis
# Perform fold-change analysis on uploaded data, unpaired
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)

# Plot fold-change analysis
mSet <- PlotFC(mSet, file.path(data_dir, "fc_0_"), "png", 72, width=NA)

# To view fold-change 
mSet$analSet$fc$fc.log


### T-Test
# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, "fdr", FALSE)

# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = file.path(data_dir, "tt_0_"), format = "png", dpi = 72, width=NA)


### Volcano Plot
# Perform the volcano analysis
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.1, TRUE, "raw")

# Create the volcano plot
mSet<-PlotVolcano(mSet,  file.path(data_dir, "volcano_0_"), 1, 0, format ="png", dpi=72, width=NA)


### ANOVA
# Perform ANOVA
mSet <- ANOVA.Anal(mSet, F, 0.05, "fisher")

# Plot ANOVA
mSet <- PlotANOVA(mSet, file.path(data_dir, "aov_0_"), "png", 72, width=NA)


### 壞掉了
# ### Correlation Analysis
# ### OPTION 1 - Heatmap specifying pearson distance and an overview
# mSet<-PlotCorrHeatMap(mSet, file.path(data_dir, "corr_0_"), "png", 72, width=NA, "col", "pearson", "bwm", "overview", F, F, 0.0)

# ### OPTION 2 - Heatmap specifying pearson correlation and a detailed view
# mSet<-PlotCorrHeatMap(mSet, file.path(data_dir, "corr_1_"),  format = "png", dpi=72, width=NA, "col", "spearman", "bwm", "detail", F, F, 999)
###




### Correlation Analysis(fixed)
source("scripts/functions/PlotCorrHeatMap_fixed.R")
 mSet <- PlotCorrHeatMap_fixed(
  mSetObj = mSet,
  out_png = "results/statistical/corr_custom.png",
  target = "col",
  method = "pearson",
  corrCutoff = 0
)


### Pattern Searching
# Perform correlation analysis on a pattern (a feature of interest in this case)
mSet<-FeatureCorrelation(mSet, "pearson", "1,6-Anhydro-beta-D-glucose")

# Plot the correlation analysis on a pattern
mSet<-PlotCorr(mSet, file.path(data_dir, "ptn_3_"), format="png", dpi=72, width=NA)


### PCA
# Perform PCA analysis
mSet<-PCA.Anal(mSet)

# Create PCA overview
mSet<-PlotPCAPairSummary(mSet, file.path(data_dir, "pca_0_"), format = "png", dpi = 72, width=NA, 5)

# Create PCA scree plot
mSet<-PlotPCAScree(mSet, file.path(data_dir, "pca_scree_0_"), "png", dpi = 72, width=NA, 5)

# Create a 2D PCA score plot
mSet<-PlotPCA2DScore(mSet, file.path(data_dir, "pca_score2d_0_"), format = "png", dpi=72, width=NA, 1, 2, 0.95, 1, 0)

# # Create a 3D PCA score plot
# mSet<-PlotPCA3DScoreImg(mSet, file.path(data_dir, "pca_score3d_0_"), "png", 72, width=NA, 1,2,3, 40)

source("scripts/functions/PlotPCA3DScoreImg_fixed.R")
mSet <- PlotPCA3DScoreImg_fixed(
  mSetObj = mSet,
  imgName = file.path(data_dir, "pca_score3d_fixed_"),
  format = "png",
  dpi = 72,
  width = NA,
  pc1 = 1,
  pc2 = 2,
  pc3 = 3,
  angle = 40
)


# Create a PCA loadings Plots
mSet<-PlotPCALoading(mSet, file.path(data_dir, "pca_loading_0_"), "png", 72, width=NA, 1,2);

# Create a PCA Biplot
mSet<-PlotPCABiplot(mSet, file.path(data_dir, "pca_biplot_0_"), format = "png", dpi = 72, width=NA, 1, 2)


### PLS-DA
mSet<-PLSR.Anal(mSet, reg=TRUE)

mSet<-PlotPLSPairSummary(mSet, file.path(data_dir, "pls_0_"), "png", 72, width=NA, 5)

mSet<-PlotPLS2DScore(mSet, file.path(data_dir, "pls_score2d_0_"), "png", 72, width=NA, 1,2,0.95,1,0)

mSet<-PlotPLS3DScoreImg(mSet, file.path(data_dir, "pls_score3d_0_"), "png", 72, width=NA, 1,2,3, 40)

mSet<-PlotPLSLoading(mSet, file.path(data_dir, "pls_loading_0_"), "png", 72, width=NA, 1, 2);

mSet<-PLSDA.CV(mSet, "5", 5,5, "Q2")


mSet<-PlotPLS.Classification(mSet, file.path(data_dir, "pls_cv_0_"), "png", 72, width=NA)

mSet<-PlotPLS.Imp(mSet, file.path(data_dir, "pls_imp_0_"), "png", 72, width=NA, "vip", "Comp. 1", 15, FALSE)

mSet<-PLSDA.Permut(mSet, 100, "accu")

mSet<-PlotPLS.Permutation(mSet, file.path(data_dir, "pls_perm_0_"), "png", 72, width=NA)
# View the 3D interactive PLS-DA score plot
mSet$imgSet$plsda.3d


### sPLS-DA
# Perform sPLS-DA analysis
mSet<-SPLSR.Anal(mSet, 5, 10, "same", "Mfold", 5, T)

# Plot sPLS-DA overview
mSet<-PlotSPLSPairSummary(mSet, file.path(data_dir, "spls_0_"), format = "png", dpi=72, width=NA, 5)

# Create 2D sPLS-DA Score Plot
mSet<-PlotSPLS2DScore(mSet, file.path(data_dir, "spls_score2d_0_"), format = "png", dpi=72, width=NA, 1, 2, 0.95, 1, 0)

# Create 3D sPLS-DA Score Plot
mSet<-PlotSPLS3DScoreImg(mSet, file.path(data_dir, "spls_score3d_0_"), format = "png", 72, width=NA, 1, 2, 3, 40)

# Create sPLS-DA loadings plot
mSet<-PlotSPLSLoading(mSet, file.path(data_dir, "spls_loading_0_"), format = "png", dpi=72, width=NA, 1,"overview")

# Perform cross-validation and plot sPLS-DA classification
mSet<-PlotSPLSDA.Classification(mSet, file.path(data_dir, "spls_cv_0_"), format = "png", dpi=72, width=NA)
# View the 3D interactive PLS-DA score plot
mSet$imgSet$splsda.3d



### orthoPLS-DA
# Perform oPLS-DA analysis
mSet<-OPLSR.Anal(mSet, reg=TRUE)

# Create a 2D oPLS-DA score plot
mSet<-PlotOPLS2DScore(mSet, file.path(data_dir, "opls_score2d_0_"), format = "png", dpi=72, width=NA, 1,2,0.95,1,0)

# Create a significant features plot
mSet<-PlotOPLS.Splot(mSet, file.path(data_dir, "opls_splot_0_"), "all", "png", 72, width=NA);

# Create a plot of features ranked by importance
mSet<-PlotOPLS.Imp(mSet, file.path(data_dir, "opls_imp_0_"), "png", 72, width=NA, "vip", "tscore", 15,FALSE)

# Create a plot of the model overview
mSet<-PlotOPLS.MDL(mSet, file.path(data_dir, "opls_mdl_0_"), format = "png", dpi=72, width=NA)

# Perform and plot oPLS-DA permutation 
mSet<-OPLSDA.Permut(mSet, 100)

mSet<-PlotOPLS.Permutation(mSet, file.path(data_dir, "opls_perm_0_"), format = "png", dpi=72, width=NA)



### SAM
# Perform SAM analysis
mSet<-SAM.Anal(mSet, "d.stat", FALSE, TRUE, 0.0, "sam_imp_0_")

# Create a SAM plot of FDR values
mSet<-PlotSAM.FDR(mSet, file.path(data_dir, "sam_fdr_0_"), format = "png", dpi=72, width=NA)

# Create a SAM plot of results
mSet<-PlotSAM.Cmpd(mSet, file.path(data_dir, "sam_imp_0_"), format = "png", dpi=72, width=NA)


### EBAM
# Perform EBAM analysis, plot EBAM analysis and create the EBAM matrix of significant features
mSet<-EBAM.Init(mSet, FALSE, TRUE, FALSE, -99.0, 0.9, file.path(data_dir, "ebam_0_"), "ebam_imp_0_")

# Create a EBAM plot of results
PlotEBAM.Cmpd(mSet, file.path(data_dir, "ebam_imp_0_"), "png", 72, width=NA)



### Hierarchical Clustering: Dendogram
# Perform hierarchical clustering and plot dendogram
mSet<-PlotHCTree(mSet, file.path(data_dir, "tree_0_"), format = "png", dpi=72, width=NA, "euclidean", "ward.D")


# ### Heatmaps
# # Perform hierarchical clustering and plot heat map
# mSet<-PlotHeatMap(mSet, file.path(data_dir, "heatmap_0_"), "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", T, T, NULL, T, F, T, T, T)

source("scripts/functions/PlotHeatMap_fixed.R")
mSet <- PlotHeatMap_fixed(
  mSetObj = mSet,
  imgName = file.path(data_dir, "heatmap_fixed.png"),
  dataOpt = "norm",
  scaleOpt = "row",
  smplDist = "euclidean",
  clstDist = "ward.D",
  palette = "bwm"
)


### K-Means
# Perform K-means analysis
mSet<-Kmeans.Anal(mSet, 3)

# Plot K-means analysis 
mSet<-PlotKmeans(mSet, file.path(data_dir, "km_0_"), format = "png", dpi=72, width=NA)


### SOM
# Perform SOM analysis
mSet<-SOM.Anal(mSet, 1, 3,"linear","gaussian")

# Plot SOM analysis
mSet<-PlotSOM(mSet, file.path(data_dir, "som_0_"), format = "png", dpi=72, width=NA)


### Random Forest
# Perform random forest analysis
mSet<-RF.Anal(mSet, 500, 7, 1)




# Plot random forest classification
mSet<-PlotRF.Classify(mSet, file.path(data_dir, "rf_0_"), format = "png", dpi=72, width=NA)

# Plot random forest variables of importance
mSet<-PlotRF.VIP(mSet, file.path(data_dir, "rf_imp_0_"), format = "png", dpi=72, width=NA)

# Plot random forest outliers 
mSet<-PlotRF.Outlier(mSet, file.path(data_dir, "rf_outlier_0_"), format = "png", dpi=72, width=NA)


### SVM
# Perform SVM 
mSet<-RSVM.Anal(mSet, 10)

mSet<-PlotRSVM.Classification(mSet, file.path(data_dir, "svm_0_"), format = "png", dpi=72, width=NA)

mSet<-PlotRSVM.Cmpd(mSet, file.path(data_dir, "svm_imp_0_"), format = "png", dpi=72, width=NA)