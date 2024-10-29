# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(viridis)  # For the inferno color palette
library(stringr)
library(readr)

# Step 1: Load normalized counts from a CSV file
counts_file <- "normalised_counts_Pirnas.csv"
normalized_counts <- read.csv(counts_file, row.names = 1)

# Step 2: Calculate CPM, log-transform, and Z-score normalize
cpm <- sweep(normalized_counts, 2, colSums(normalized_counts), FUN = "/") * 1e6
logCPM <- log2(cpm + 1)  # Log-transform with pseudo-count of 1
logCPM_z <- t(scale(t(logCPM)))  # Z-score normalization across samples for each gene

# List of specific row names (genes) to match
row_names_to_match <- c(
  "dre-miR-184",
  "dre-miR-10b-5p",
  "dre-miR-129-5p",
  "dre-miR-183-5p",
  "dre-miR-181a-5-3p",
  "dre-miR-193b-5p",
  "dre-miR-200a-5p")

# Function to match row names
match_row <- function(row_name) {
  matched_row <- rownames(logCPM_z)[grep(row_name, rownames(logCPM_z))]
  return(matched_row)
}

# Match each row name in the list and extract the data
matched_rows <- sapply(row_names_to_match, match_row, simplify = FALSE)
matched_rows_vector <- unlist(matched_rows)
heatmap_data <- logCPM_z[matched_rows_vector, ]

# Identify missing rows (for debugging, if needed)
missing_rows <- setdiff(row_names_to_match, matched_rows_vector)
if (length(missing_rows) > 0) {
  cat("The following rows are missing from logCPM_z and won't appear in the heatmap:\n")
  print(missing_rows)
}

# Modify row names to the format CHROM:<chromosome>:<start-end>
new_row_names <- rownames(heatmap_data)
for (i in seq_along(new_row_names)) {
  # Match the chromosome number and start-end positions
  if (grepl(".*?:(\\d+):.*?_(\\d+-\\d+)", new_row_names[i])) {
    new_row_names[i] <- sub(".*?:(\\d+):.*?_(\\d+-\\d+).*", "CHROM:\\1:\\2", new_row_names[i])
  } else if (grepl("^[0-9-]+$", new_row_names[i])) {
    # If numeric-only, prefix with CHROM (for simpler row names without complex structure)
    new_row_names[i] <- paste("CHROM:", new_row_names[i], sep = "")
  }
}
rownames(heatmap_data) <- new_row_names

# Remove unwanted duplicates (if any)
rows_to_remove <- c(
  "none"
)
heatmap_data_filtered <- heatmap_data[!(rownames(heatmap_data) %in% rows_to_remove), ]

# Define custom column labels for the heatmap
custom_column_labels <- c("High.2_1", "High.2_2", "High.2_3", "High.2_4", "High.2_5", "High.2_6", "High.2_7", "High.2_8", "High.2_9", "High.2_20",
                          "Low.1_1", "Low.1_2", "Low.1_3", "Low.1_4", "Low.1_5", "Low.1_6","Low.1_7", "Low.1_8", "Low.1_9", "Low.1_10",
                          "Low.2_1", "Low.2_2", "Low.2_3", "Low.2_4", "Low.2_5", "Low.2_6","Low.2_7", "Low.2_8", "Low.2_9", "Low.2_10",
                          "High.1_1", "High.1_2", "High.1_3", "High.1_4", "High.1_5", "High.1_6", "High.1_7", "High.1_8", "High.1_9", "High.1_10")

# Verify that the column names in heatmap_data_filtered match custom_column_labels
colnames(heatmap_data_filtered) <- custom_column_labels  # Ensure columns have custom names

# Create treatment labels for annotation, matching custom_column_labels order
treatment_labels <- c("High2", "High2", "High2", "High2", "High2", "High2", "High2", "High2", "High2", "High2",
                      "Low1", "Low1", "Low1", "Low1", "Low1", "Low1","Low1", "Low1", "Low1", "Low1",
                      "Low2", "Low2", "Low2", "Low2", "Low2", "Low2", "Low2", "Low2", "Low2", "Low2",
                      "High1", "High1", "High1", "High1", "High1", "High1", "High1", "High1", "High1", "High1")

# Define the annotation data frame, ensuring the row names match custom_column_labels
annotation_col <- data.frame(Treatment = factor(treatment_labels, levels = c("High1", "High2", "Low1", "Low2")))
rownames(annotation_col) <- custom_column_labels

# Define colors for each treatment group in annotation_colors
annotation_colors <- list(
  Treatment = c("High1" = "#E69F00", "High2" = "#56B4E9", "Low1" = "#009E73", "Low2" = "#F0E442")
)

# Plot the heatmap with filtered data and annotations, without clustering columns
pheatmap(heatmap_data_filtered, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Do not cluster columns; keep original order
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontface = "bold",
         scale = "none",  # Data is already z-scored
         color = inferno(256),
         labels_col = custom_column_labels, 
         annotation_col = annotation_col,  # Add the column annotations for treatments
         annotation_colors = annotation_colors,  # Apply colors for each treatment group
         annotation_legend = TRUE,  # Show legend for the annotations
         main = "Z-Score Normalized logCPM for Selected miRNAs with Treatment Grouping")
