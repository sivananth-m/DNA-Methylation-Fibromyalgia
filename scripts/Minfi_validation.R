getwd()
setwd("D:/Endo")
# ==============================================================================
# MINFI WORKFLOW FOR FIBROMYALGIA METHYLATION ANALYSIS
# ==============================================================================

# 1. Setup Environment and Load Libraries ----------------------------------------

# Install packages if needed (uncomment and run these lines once)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("minfi", "DMRcate", "limma", "RColorBrewer", "missMethyl"))

# Load necessary libraries
library(minfi)        # Core package for 450k/EPIC array analysis
library(limma)        # For linear modeling and differential analysis
library(DMRcate)      # For finding differentially methylated regions (optional but recommended)
library(RColorBrewer) # For plotting colors
library(missMethyl)   # For functional normalization and QC
library(doParallel)   # For parallel processing (speeds things up)
library(ggplot2)      # For enhanced visualizations
packageVersion("minfi")


# Set up parallel processing (adjust cores as needed)
registerDoParallel(cores = detectCores() - 1)
message(paste("Using", getDoParWorkers(), "cores for parallel processing."))


# 2. Data Loading and Phenotype Acquisition --------------------------------------

# NOTE: You MUST place your raw IDAT files (e.g., GSM2396835_R01C01_Grn.idat)
# into a folder named 'idat_files' in the same directory as this script.

# Define the directory where your IDAT files are stored
baseDir <- "D:/Endo"

# Read the sample sheet (now using the specified SampleSheet.csv file)
targets <- read.metharray.sheet(baseDir, pattern = "SampleSheet.csv")

# --- START OF CRITICAL UPDATES FOR FIBROMYALGIA DATA ---
# Check if the 'Sample_Group' column exists (based on your new file)
if (!"Sample_Group" %in% names(targets)) {
  stop("Error: The column 'Sample_Group' was not found in the Sample Sheet. Please check the column headers in your SampleSheet.csv file.")
}

# minfi functions often expect a column named 'Group' for phenotype grouping
targets$Group <- targets$Sample_Group

# Ensure the 'Group' column is a factor, using the correct levels: "Control" and "Fibromyalgia".
targets$Group <- factor(targets$Group, levels = c("Pre_Receptive", "Receptive"))
# --- END OF CRITICAL UPDATES ---

print("--- Targets (Sample Sheet) Loaded ---")
print(targets)


# 3. Load Raw Data (RGChannelSet) ------------------------------------------------

# Load the raw intensity data from the IDAT files specified in the targets object.
# This creates an RGChannelSet, which holds the Red and Green raw intensity channels.
rgSet <- read.metharray.exp(targets = targets, force = TRUE)

# Add sample names to the RGSet (optional, but helpful for plotting)
sampleNames(rgSet) <- targets$Sample_Name
print("--- Raw Data (RGChannelSet) Loaded ---")


# 4. Quality Control (QC) --------------------------------------------------------

# Plot the distribution of intensity values before normalization.
# NOTE: This is the correct, two-step approach for QC plotting:
# 1. qc(rgSet) calculates the quality control metrics.
# 2. plotQC() plots these metrics (Median intensity vs. Control Probes).
# 4. Quality Control (QC)
# ---- QC on raw data (before normalization) ----
qcReport(rgSet, pdf = TRUE)
print("QC Report (pre-normalization) generated as qcReport.pdf")

# ---- Functional normalization (keep this as-is) ----
mSetSq <- preprocessFunnorm(rgSet)
print("--- Data Preprocessed using preprocessFunnorm ---")

# Extract the Methylated and Unmethylated intensities before conversion
mSetRaw <- preprocessRaw(rgSet)

qcData <- getQC(mSetRaw)
plotQC(qcData)

# Calculate detection p-values (to identify probes with low confidence)
detP <- detectionP(rgSet)
print(paste("Detection P-value matrix dimension:", dim(detP)[1], "probes,", dim(detP)[2], "samples."))


# 5. Normalization and Preprocessing ---------------------------------------------

# Convert RGChannelSet to MethylSet, then to GenomicRatioSet (for functional normalization).
# We use preprocessFunnorm, which is highly recommended for studies with
# cell composition heterogeneity or differences in tissue/group, like this one.
mSetSq <- preprocessFunnorm(rgSet)
print("--- Data Preprocessed using preprocessFunnorm ---")

# Apply QC filtering to the normalized data:
# a) Filter probes with high detection p-values (>0.01) in any sample
#    (Note: This is a strict filter; adjust threshold if needed)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSq.filtered <- mSetSq[keep, ]

# b) Filter out probes known to be cross-reactive or containing common SNPs
#    The minfi/missMethyl packages provide built-in lists for this.
#    Get the relevant annotation (e.g., 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19')
#    Replace '450k' with 'EPIC' if you are using the EPIC array.
if (rgSet@annotation["array"] == "450k") {
  ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  # Filter out probes on sex chromosomes (X/Y)
  keepCpGs <- !(featureNames(mSetSq.filtered) %in% ann$Name[ann$chr %in% c("chrX", "chrY")])
  mSetSq.filtered <- mSetSq.filtered[keepCpGs, ]
  
  # Filter out cross-reactive probes and SNP-associated probes
  data(snps.450k.hg19)
  mSetSq.filtered <- dropLociWithSnps(mSetSq.filtered, snps = snps.450k.hg19)
  
} else if (rgSet@annotation["array"] == "EPIC") {
  ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  # Filter out probes on sex chromosomes (X/Y)
  keepCpGs <- !(featureNames(mSetSq.filtered) %in% ann$Name[ann$chr %in% c("chrX", "chrY")])
  mSetSq.filtered <- mSetSq.filtered[keepCpGs, ]
  
  # Filter out cross-reactive probes and SNP-associated probes
  data(snps.epic.hg19)
  mSetSq.filtered <- dropLociWithSnps(mSetSq.filtered, snps = snps.epic.hg19)
}

# Check remaining probes
print(paste("Probes remaining after filtering:", nrow(mSetSq.filtered)))


# 6. Differential Methylation Analysis (DMP) --------------------------------------

# Convert to Beta (methylation) and M-values (better for statistical analysis)
beta.values <- getBeta(mSetSq.filtered)
m.values <- getM(mSetSq.filtered)
pheno <- pData(mSetSq.filtered)

# A. Using minfi's dedicated function (DMPfind)
# This is a fast way to find single-CpG DMPs using an empirical Bayes method.
dmp <- dmpFinder(m.values, pheno = pheno$Group, type = "categorical")

# Calculate the mean difference in Beta values: Fibromyalgia - Control
dmp$abs_diff_beta <- rowMeans(beta.values[, pheno$Group == "Receptive"]) - rowMeans(beta.values[, pheno$Group == "Pre_Receptive"])

# Sort by p-value
dmp.sorted <- dmp[order(dmp$pval), ]

# Display top 10 differentially methylated positions (DMPs)
print("--- Top 10 Differentially Methylated Positions (DMPs) ---")
print(head(dmp.sorted, 10))
top10_DMPs<-(head(dmp.sorted, 10))
write.csv(top10_DMPs, file = "Top10_DMPs.csv", row.names = TRUE)



# B. Using the standard limma approach (recommended for more complex designs)
# Create the design matrix
design <- model.matrix(~Group, data = pheno)

# Fit the linear model to the M-values
fit <- lmFit(m.values, design)

# Perform empirical Bayes smoothing (essential for good results)
fit.bayes <- eBayes(fit)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


# Extract differential methylation results for the 'Fibromyalgia' vs 'Control' contrast
# The coefficient (column 2) corresponds to the 'GroupFibromyalgia' comparison
DMP.limma <- topTable(fit.bayes, coef = 2, number = Inf, genelist = ann[rownames(fit$coefficients), ])
DMP.limma$abs_diff_beta <- dmp.sorted[rownames(DMP.limma), "abs_diff_beta"]

print("--- Top 10 DMPs using limma (M-values) ---")
print(head(DMP.limma, 10))

# Save results
write.csv(DMP.limma, file = "minfi_DMP_results_limma.csv", row.names = TRUE)
print("Results saved to 'minfi_DMP_results_limma.csv'.")


# 7. Visualization ---------------------------------------------------------------

# A. MDS Plot (Multi-Dimensional Scaling)
# This checks for clustering based on phenotype/batch effects after normalization.
betaVals <- getBeta(mSetSq.filtered)
pca <- prcomp(t(betaVals), scale. = TRUE)

# Plot first 2 PCs
plot(pca$x[,1], pca$x[,2],
     col = brewer.pal(8, "Dark2")[factor(pheno$Group)],
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA Plot (Post-Normalization)")
legend("topright", legend = levels(factor(pheno$Group)),
       col = brewer.pal(8, "Dark2"), pch = 19)

# B. Density Plot of Beta Values
# Should show a bimodal distribution (unmethylated ~0, methylated ~1)
par(mfrow=c(1,1))
densityPlot(beta.values, main = "Beta Value Distribution (Post-Normalization)",
            pal = brewer.pal(8, "Dark2")[factor(pheno$Group)], legend = FALSE)
legend("top", legend = levels(factor(pheno$Group)),
       text.col = brewer.pal(8, "Dark2")[1:length(unique(pheno$Group))],
       bty = "n", cex = 1)


# C. Volcano Plot (to visualize DMPs)
# Volcano Plot for limma results (M-values)
DMP.limma$logFC <- DMP.limma$logFC / 0.59 # Scale M-value logFC to roughly reflect Beta differences
DMP.limma$is_sig <- ifelse(DMP.limma$adj.P.Val < 0.05 & abs(DMP.limma$logFC) > 0.5, "Significant", "Non-Significant")

p_volcano <- ggplot(DMP.limma, aes(x = logFC, y = -log10(adj.P.Val), color = is_sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "#E53E3E", "Non-Significant" = "#9EAEB9")) +
  labs(x = "Log Fold Change (M-Value Equivalent)",
       y = "-log10(Adjusted P-value)",
       title = "Volcano Plot of Differential Methylation in Fibromyalgia") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())

print(p_volcano)


# 8. Differential Methylated Region (DMR) Analysis (Optional) ---------------------

# DMR analysis looks at contiguous probes to find regions, which is often more biologically
# meaningful than single-CpG analysis.

# Calculate the correlation matrix (optional, but improves DMRcate results)
# The annotation object 'ann' needs to be defined from section 5.
library(missMethyl)
library(DMRcate)

# Step 1: Annotate CpGs
cpg.groups <- cpg.annotate(
  datatype = "array",
  object = mSetSq.filtered,
  arraytype = "450K",
  analysis.type = "differential",
  design = design,
  coef = 2
)

# Step 2: Run DMRcate
dmrcate.results <- dmrcate(cpg.groups, pcut = 0.05, lambda = 1000, C = 2)

# Convert DMResults object to GRanges
dmr.ranges <- extractRanges(dmrcate.results)

# Convert to a dataframe
dmr.table <- as.data.frame(dmr.ranges)

# Save to CSV
write.csv(dmr.table, "DMRcate_results.csv", row.names = FALSE)




# ==============================================================================
# END OF SCRIPT
# ==============================================================================
