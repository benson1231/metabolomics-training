library(MetaboAnalystR)

# Initialize Objects
mSet <- InitDataObjects("nmrpeak", "stat", FALSE);

# Fetch data and unzipping
download.file("https://rest.xialab.ca/api/download/metaboanalyst/nmr_peaks.zip", destfile = "nmr_peaks.zip", method = "auto");
UnzipUploadedFile("nmr_peaks.zip", "upload");

# Read data in 
mSet <- Read.PeakList(mSet, "upload");
  
# Perform grouping of peaks
mSet <- GroupPeakList(mSet, 0.025, 30.0);
  
# Form peak groups
mSet <- SetPeakList.GroupValues(mSet)

# Run the sanity check, it will return a series of messages if the data is suitable for subsequent analyses. 
mSet <- SanityCheckData(mSet)

# Replace missing/zero values with a minimum positive value
mSet <- ReplaceMin(mSet)

# View messages collected during ReplaceMin()
mSet$msgSet$replace.msg

######### Alternative Step 2: Replace missing values with KNN imputed values   
mSet <- ImputeMissingVar(mSet, method="knn_smp")

# Check if the sample size is too small, returns a 0 if the data passes the check
mSet<-IsSmallSmplSize(mSet)

### OPTION 1) Perform Probabilistic Quotient Normalization based upon a reference sample
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "SamplePQN", "NULL", "NULL", "P037", ratio=FALSE, ratioNum=20)

# View feature normalization
mSet<-PlotNormSummary(mSet, "feature_norm", format="png", dpi=72, width=NA)

# View sample normalization
mSet<-PlotSampleNormSummary(mSet, "sample_norm", format="pdf", width=NA)