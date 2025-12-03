# ==============================================================================
# Script: 02_Plotting.R
# Purpose: Generate visualization for Fibromyalgia Methylation Analysis
# Author: [Your Name]
# ==============================================================================

# 1. Load Required Libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
library(ggplot2)
library(dplyr)

# ==============================================================================
# Part 1: Load Results Data
# ==============================================================================
# We load the results file we saved in the previous step
dmp_file <- "DMP_Results_Fibromyalgia_vs_Control.csv"

if(file.exists(dmp_file)){
    results <- read.csv(dmp_file)
} else {
    stop("DMP Results file not found! Check your working directory.")
}

# ==============================================================================
# Part 2: Volcano Plot
# ==============================================================================
print("Generating Volcano Plot...")

# Add a column to categorize significance
# Significant = Adj.P.Val < 0.05
# Hyper-methylated (Up) = deltaBeta > 0
# Hypo-methylated (Down) = deltaBeta < 0
results$Category <- "Not Significant"
results$Category[results$adj.P.Val < 0.05 & results$deltaBeta > 0] <- "Hyper-methylated"
results$Category[results$adj.P.Val < 0.05 & results$deltaBeta < 0] <- "Hypo-methylated"

# Create the Plot
volcano_plot <- ggplot(results, aes(x = deltaBeta, y = -log10(P.Value), color = Category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Hyper-methylated" = "#E41A1C", 
                                "Hypo-methylated" = "#377EB8", 
                                "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.0001), linetype = "dashed", color = "black") + # Visual threshold line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot: Differential Methylation",
       subtitle = "Fibromyalgia vs. Control",
       x = "Delta Beta (Methylation Difference)",
       y = "-log10(P-Value)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save the plot
ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 6)


# ==============================================================================
# Part 3: Manhattan Plot
# ==============================================================================
print("Generating Manhattan Plot...")

# Prepare data: Convert numeric Chromosomes and calculate cumulative position
# Filter out X/Y for cleaner plot if needed, or keep them.
manhattan_data <- results %>%
    filter(CHR %in% c(1:22)) %>%  # focusing on autosomes for speed
    mutate(CHR = as.numeric(CHR)) %>%
    arrange(CHR, MAPINFO) %>%
    mutate(BPcum = NA)

# Calculate cumulative position for plotting continuous x-axis
data_cum <- manhattan_data %>%
    group_by(CHR) %>%
    summarise(max_bp = max(MAPINFO)) %>%
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
    select(CHR, bp_add)

manhattan_data <- manhattan_data %>%
    inner_join(data_cum, by = "CHR") %>%
    mutate(bp_cum = MAPINFO + bp_add)

# Calculate axis labels (center of each chromosome)
axis_set <- manhattan_data %>%
    group_by(CHR) %>%
    summarize(center = mean(bp_cum))

# Create the Plot
manhattan_plot <- ggplot(manhattan_data, aes(x = bp_cum, y = -log10(P.Value))) +
    # Show all points in grey
    geom_point(aes(color = as.factor(CHR)), alpha = 0.5, size = 1.0) +
    scale_color_manual(values = rep(c("#2c3e50", "#e74c3c"), 22)) + # Alternating colors
    # Custom X axis
    scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(title = "Manhattan Plot: Genome-wide Significant Probes",
         x = "Chromosome",
         y = "-log10(P-Value)")

# Save the plot
ggsave("manhattan_plot.png", plot = manhattan_plot, width = 12, height = 6)


# ==============================================================================
# Part 4: TNXB Gene Boxplot
# ==============================================================================
print("Generating TNXB Boxplot...")

# NOTE: This section requires 'myNorm' and 'myLoad' objects from the ChAMP pipeline.
# If they are missing, this part will be skipped to prevent errors.

if(exists("myNorm") && exists("myLoad")) {
    
    # 1. Find the most significant probe for TNXB
    # We grep for TNXB in the results we loaded
    tnxb_probes <- results[grep("TNXB", results$gene), ]
    
    # Check if we found any probes
    if(nrow(tnxb_probes) > 0) {
        # Sort by P-value to get the best one
        tnxb_probes <- tnxb_probes %>% arrange(P.Value)
        top_probe_id <- tnxb_probes$X[1] # Assuming first column is Probe ID (or check rowname)
        
        # If ID is not in column X, try to use row names if results was a matrix converted
        if(is.null(top_probe_id)) top_probe_id <- rownames(tnxb_probes)[1]
        
        print(paste("Plotting probe:", top_probe_id))
        
        # 2. Extract Beta Values for this probe
        # Corresponds to: myNorm[probe_id, ]
        probe_data <- data.frame(
            Beta = myNorm[top_probe_id, ],
            Group = myLoad$pd$Sample_Group
        )
        
        # 3. Create Boxplot
        tnxb_plot <- ggplot(probe_data, aes(x = Group, y = Beta, fill = Group)) +
            geom_boxplot(outlier.shape = NA, alpha = 0.7) +
            geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
            scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
            theme_minimal() +
            labs(title = paste("Methylation of TNXB (", top_probe_id, ")"),
                 subtitle = "Linked to Ehlers-Danlos / Connective Tissue Integrity",
                 y = "Beta Value (Methylation Level)",
                 x = "Group") +
            theme(legend.position = "none",
                  plot.title = element_text(face = "bold", hjust = 0.5))
        
        # Save
        ggsave("TNXB_Boxplot.png", plot = tnxb_plot, width = 6, height = 6)
        print("TNXB Boxplot saved.")
        
    } else {
        print("No TNXB probes found in the results file.")
    }
    
} else {
    print("WARNING: 'myNorm' or 'myLoad' objects not found. Skipping Boxplot.")
    print("To generate the boxplot, you must run the ChAMP Load and Norm steps first.")
}

print("All plots generated successfully!")
