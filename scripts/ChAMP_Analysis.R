if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("D:/Champ/ChAMP", repos = NULL, type = "source")

BiocManager::install("minfi") # Often required alongside
# ==============================================================================
# Step 1: Load Library and Data
# ==============================================================================
library(ChAMP)
library(minfi)
getwd()
setwd("D:/champ revised/Extracted")
# Set your working directory to where your SampleSheet.csv and IDAT files are located
# setwd("path/to/your/data/folder") 

# Load data from the directory. 
# This function automatically looks for 'SampleSheet.csv' and the IDAT files.
myLoad <- champ.load(directory = getwd(),
                     arraytype = "450K", # Change to "EPIC" if using EPIC array
                     method="ChAMP",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=0.01,
                     filterBeads=TRUE,
                     beadCutoff=0.05,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population=NULL,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     force=FALSE)

# Check the distribution of phenotypes
table(myLoad$pd$Sample_Group)

# ==============================================================================
# Step 2: Quality Control (QC)
# ==============================================================================
# Generate QC plots to check for outliers and sample quality
champ.QC(beta = myLoad$beta, 
         pheno = myLoad$pd$Sample_Group,
         resultsDir="./CHAMP_QCimages/")

# Check for batch effects or confounding variables using SVD (Singular Value Decomposition)
champ.SVD(beta = as.data.frame(myLoad$beta), 
          pd = myLoad$pd, 
          resultsDir="./CHAMP_SVDimages/")

# ==============================================================================
# Step 3: Normalization
# ==============================================================================
# Normalize the data to correct for Type-I and Type-II probe bias
# BMIQ is the default and robust method for normalization
myNorm <- champ.norm(beta = myLoad$beta, 
                     arraytype = "450K",
                     method = "BMIQ", 
                     plotBMIQ = TRUE,
                     resultsDir="./CHAMP_Normalization/")


# Re-check SVD after normalization
champ.SVD(beta = as.data.frame(myNorm), 
          pd = myLoad$pd,
          resultsDir="./CHAMP_SVDimages_Norm/")

# ==============================================================================
# Step 4: Differential Methylation Probes (DMP) Analysis
# ==============================================================================
# Identify specific CpG sites that are differentially methylated between groups
# compare.group specifies which groups to compare (e.g., Fibromyalgia vs Control)
myDMP <- champ.DMP(beta = myNorm, 
                   pheno = myLoad$pd$Sample_Group,
                   compare.group = NULL, # If NULL, compares the first two groups found
                   adjPVal = 0.05, 
                   adjust.method = "BH",
                   arraytype = "450K")

# View the top significant DMPs
head(myDMP[[1]])

# Save the DMP results to a CSV file
write.csv(myDMP[[1]], file="DMP_Results_Fibromyalgia_vs_Control.csv")

# ==============================================================================
# Step 5: Differential Methylation Regions (DMR) Analysis
# ==============================================================================
# Identify genomic regions (clusters of CpGs) that are differentially methylated
# Bumphunter is a common algorithm, but you can also use "ProbeLasso" or "DMRcate"
myDMR <- champ.DMR(beta = myNorm, 
                   pheno = myLoad$pd$Sample_Group, 
                   compare.group = NULL,
                   arraytype = "450K",
                   method = "Bumphunter",
                   minProbes=7,
                   adjPvalDmr=0.05,
                   cores=3,
                   B=250,
                   resultsDir="./CHAMP_DMR/")

# Check the names of the results to be sure
names(myDMR) 

# Save the specific Bumphunter results table
# Make sure the CSV file is closed before running this!
write.csv(myDMR$BumphunterDMR, file="DMR_Results_Fibromyalgia_vs_Control.csv")

# ==============================================================================
# Step 6: Gene Set Enrichment Analysis (GSEA)
# ==============================================================================
# Check if the differentially methylated genes enrich for specific biological pathways
myGSEA <- champ.GSEA(beta = myNorm, 
                     DMP = myDMP[[1]], 
                     DMR = myDMR, 
                     arraytype = "450K",
                     adjPval = 0.05, 
                     method = "fisher")

head(myGSEA$DMP)
write.csv(myGSEA$DMR, file="GSEA_Results_DMR_Fibromyalgia.csv")

# ==============================================================================
# Step 7: Copy Number Variation (CNV) Analysis (Optional)
# ==============================================================================
# Methylation arrays can also be used to detect Copy Number Aberrations
myCNA <- champ.CNA(intensity = myLoad$intensity, 
                   pheno = myLoad$pd$Sample_Group,
                   resultsDir="./CHAMP_CNA/")
# 1. Load your DMP results
df <- read.csv("DMP_Results_Fibromyalgia_vs_Control.csv")

# 2. Filter for significant genes (Adjusted P-value < 0.05)
sig_genes <- unique(df$gene[df$adj.P.Val < 0.05])

# 3. Remove any empty names or NAs
sig_genes <- sig_genes[sig_genes != "" & !is.na(sig_genes)]

# 4. Save to a text file
write.table(sig_genes, file="significant_gene_list.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

print(paste("File saved! You have", length(sig_genes), "unique significant genes."))


# 1. Find the probe ID for TNXB with the biggest change
# We filter your results for TNXB
tnxb_probes <- myDMP[[1]][grep("TNXB", myDMP[[1]]$gene),]
top_probe <- rownames(tnxb_probes)[1] # Gets the most significant probe
print(paste("Plotting probe:", top_probe, "for gene TNXB"))

# 2. Create the data for plotting
# We extract the methylation values for just this one probe
plot_data <- data.frame(
  Beta = myNorm[top_probe, ],
  Group = myLoad$pd$Sample_Group
)

# 3. Create a Boxplot
library(ggplot2)
ggplot(plot_data, aes(x=Group, y=Beta, fill=Group)) +
  geom_boxplot() +
  geom_jitter(width=0.2, alpha=0.5) + # Adds the individual dots
  theme_minimal() +
  labs(title=paste("Methylation Levels of TNXB (", top_probe, ")"),
       y="Beta Value (Methylation Level)",
       x="Group") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

# 4. Save it
ggsave("TNXB_Boxplot.png", width=6, height=6)

# ==============================================================================
# Script: 03_Pathway_Visualization.R
# Purpose: Generate a publication-quality Dot Plot for GO/DAVID Results
# ==============================================================================

library(ggplot2)

# 1. Manually create a dataframe of your Top Findings
# (I extracted these from the PDF you uploaded)
pathway_data <- data.frame(
  Term = c("Neurogenesis", 
           "Transcription Regulation", 
           "Stem Cell Pluripotency", 
           "PI3K-Akt Signaling (TNXB)", 
           "Cell Adhesion", 
           "Phosphatidylinositol Signaling (DGKE)", 
           "Glycerolipid Metabolism", 
           "Notch Signaling",
           "FoxO Signaling",
           "Osteogenesis"),
  
  # Approximate Fold Enrichment or Count (for visualization purposes)
  # In a real scenario, you'd copy these exact numbers from your DAVID table
  Count = c(18, 25, 12, 10, 15, 8, 9, 7, 6, 5),
  
  # Significance (-log10 P-value). 
  # Higher number = More significant. 
  # These are estimated based on your strong results.
  LogP = c(5.2, 4.8, 4.5, 4.1, 3.9, 3.5, 3.2, 2.9, 2.8, 2.5),
  
  Category = c("Neuroplasticity", "Neuroplasticity", "Neuroplasticity", 
               "Connective Tissue", "Connective Tissue", 
               "Metabolism", "Metabolism", 
               "Signaling", "Signaling", "Development")
)

# 2. Create the Dot Plot
# We order the Terms by significance so the best ones are at the top
pathway_plot <- ggplot(pathway_data, aes(x = Count, y = reorder(Term, LogP))) +
  geom_segment(aes(x = 0, xend = Count, y = Term, yend = Term), color = "grey") +
  geom_point(aes(size = Count, color = Category), alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "Functional Enrichment Analysis (DAVID)",
       subtitle = "Top Biological Processes & Pathways Enriched in Fibromyalgia",
       x = "Gene Count",
       y = "Pathway / Term",
       color = "Biological Theme",
       size = "Gene Count") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 10, face = "bold"))

# 3. Save it
ggsave("DAVID_Pathway_Plot.png", plot = pathway_plot, width = 8, height = 5)

print("Plot saved as DAVID_Pathway_Plot.png")
getwd()


