# ==============================================================================
# Script: Consensus Gene Discovery (Minfi vs ChAMP)
# ==============================================================================

# 1. Load Files
minfi_dmrs <- read.csv("filtered_DMR.csv")
champ_dmps <- read.csv("DMP_Results_Fibromyalgia_vs_Control.csv")

# 2. Extract Genes from Minfi (Handling comma-separated lists)
# Split genes if multiple are listed (e.g., "GeneA,GeneB")
minfi_gene_list <- unlist(strsplit(as.character(minfi_dmrs$overlapping.genes), ","))
minfi_gene_list <- unique(trimws(minfi_gene_list)) # Remove spaces/duplicates

# 3. Extract Genes from ChAMP (Significant only)
champ_sig <- champ_dmps[champ_dmps$adj.P.Val < 0.05, ]
champ_gene_list <- unique(champ_sig$gene)

# 4. Find the "Gold Standard" Overlap
consensus_genes <- intersect(minfi_gene_list, champ_gene_list)

# 5. Save the list
write.table(consensus_genes, "Consensus_Genes_Minfi_ChAMP.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

print(paste("Number of Consensus Genes:", length(consensus_genes)))
print(consensus_genes)
library(VennDiagram)

# Create the Venn
venn.plot <- draw.pairwise.venn(
  area1 = length(minfi_gene_list),
  area2 = length(champ_gene_list),
  cross.area = length(consensus_genes),
  category = c("Minfi (Genes)", "ChAMP (Genes)"),
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  cat.pos = c(0, 0),
  scaled = FALSE # Important because ChAMP is huge compared to Minfi
)

# Save
png("results/Gene_Overlap_Venn.png")
grid.draw(venn.plot)
dev.off()