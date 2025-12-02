
# ==============================================================================
# Project: Epigenetic Signatures of Fibromyalgia
# Pipeline: ChAMP (Illumina 450k/EPIC)
# Author: Sivananth M.
# ==============================================================================

# 1. Load Libraries
library(ChAMP)
library(minfi)

# 2. Load Data
# Assumes SampleSheet.csv and IDAT files are in the working directory
myLoad <- champ.load(directory = getwd(),
                     arraytype = "450K",
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

# 3. Quality Control
champ.QC(beta = myLoad$beta, 
         pheno = myLoad$pd$Sample_Group,
         resultsDir="./CHAMP_QCimages/")

# 4. Normalization (BMIQ)
myNorm <- champ.norm(beta = as.data.frame(myLoad$beta), 
                     arraytype = "450K",
                     method = "BMIQ", 
                     plotBMIQ = TRUE,
                     resultsDir="./CHAMP_Normalization/")

# 5. SVD Analysis (Batch Effects)
champ.SVD(beta = as.data.frame(myNorm), 
          pd = myLoad$pd,
          resultsDir="./CHAMP_SVDimages_Norm/")

# 6. Differential Methylation Probes (DMP)
myDMP <- champ.DMP(beta = myNorm, 
                   pheno = myLoad$pd$Sample_Group,
                   compare.group = NULL,
                   adjPVal = 0.05, 
                   adjust.method = "BH",
                   arraytype = "450K")

# Save DMP Results
write.csv(myDMP[[1]], file="DMP_Results_Fibromyalgia_vs_Control.csv")

# 7. Differential Methylation Regions (DMR)
# Using Bumphunter method
myDMR <- champ.DMR(beta = as.data.frame(myNorm), 
                   pheno = myLoad$pd$Sample_Group, 
                   compare.group = NULL,
                   arraytype = "450K",
                   method = "Bumphunter",
                   minProbes = 7,
                   adjPvalDmr = 0.05,
                   cores = 3,
                   B = 250,
                   resultsDir = "./CHAMP_DMR/")

# Save DMR Results
write.csv(myDMR$BumphunterDMR, file="DMR_Results_Fibromyalgia_vs_Control.csv")

# 8. Gene Set Enrichment Analysis (GSEA)
myGSEA <- champ.GSEA(beta = as.data.frame(myNorm), 
                     DMR = myDMR, 
                     arraytype = "450K",
                     adjPval = 0.05, 
                     method = "fisher")

# Save GSEA Results
write.csv(myGSEA$DMR, file="GSEA_Results_DMR_Fibromyalgia.csv")
