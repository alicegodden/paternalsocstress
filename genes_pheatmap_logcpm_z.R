# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(viridis)  # For the inferno color palette
library(stringr)
library(readr)

# Step 1: Load normalized counts from a CSV file
counts_file <- "willian_rna_norm_counts_deseq2.csv"
normalized_counts <- read.csv(counts_file, row.names = 1)

# Step 2: Calculate CPM, log-transform, and Z-score normalize
cpm <- sweep(normalized_counts, 2, colSums(normalized_counts), FUN = "/") * 1e6
logCPM <- log2(cpm + 1)  # Log-transform with pseudo-count of 1
logCPM_z <- t(scale(t(logCPM)))  # Z-score normalization across samples for each gene


# Load the CSV file
data <- read.csv("mir_matched_counts_10b.csv")

# Replace `Gene` with the actual column name containing your gene names
row_names_to_match <- data$gene


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
custom_column_labels <- c("E100.1", "E100.2", "E100.3", "E100.4", "E101.1",
                          "E101.2", "E101.3", "E101.4", "E102.1", "E102.2",
                          "E102.3", "E102.4", "E103.1", "E103.2", "E103.3", 
                          "E103.4", "E104.1", "E104.2", "E104.3", "E104.4", 
                          "E105.1", "E105.2", "E105.3", "E105.4", "E106.1", 
                          "E106.2", "E106.3", "E106.4", "E107.1", "E107.2",
                          "E107.3", "E107.4", "E108.1", "E108.2", "E108.3",
                          "E108.4", "E109.1", "E109.2", "E109.3", "E109.4", 
                          "E110.1", "E110.2", "E110.3", "E110.4", "E111.1", 
                          "E111.2", "E111.3", "E111.4", "E112.4", "E112.5", 
                          "E112.7", "E112.8", "E113.10", "E113.11", "E113.13", 
                          "E113.9", "E114.1", "E114.3", "E114.4", "E114.5",
                          "E115.1", "E115.2", "E115.4", "E115.5", "E80.1",
                          "E80.2", "E80.3", "E80.4", "E81.1", "E81.2", 
                          "E81.4", "E81.5", "E82.1", "E82.2", "E82.3", 
                          "E82.5", "E83.1", "E83.3", "E83.4", "E83.5")

# Verify that the column names in heatmap_data_filtered match custom_column_labels
colnames(heatmap_data_filtered) <- custom_column_labels  # Ensure columns have custom names

# Create treatment labels for annotation, matching custom_column_labels order
treatment_labels <- c("High", "High", "High", "High", "High", "High", "High", "High",
                      "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", 
                      "High", "High", "High", "High", "High", "High", "High", "High",
                      "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", 
                      "Low", "Low", "Low", "Low", "Low", "Low", "Low", "High", "High",
                      "High", "High", "Low", "Low", "Low", "Low", "High", "High", "High",
                      "High", "High", "High", "High", "High", "Low", "Low", "Low", "Low",
                      "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low",
                      "Low", "Low", "High", "High", "High", "High", "High", "High",
                      "High", "High")

# Define the annotation data frame, ensuring the row names match custom_column_labels
annotation_col <- data.frame(Treatment = factor(treatment_labels, levels = c("High", "Low")))
rownames(annotation_col) <- custom_column_labels

# Define colors for each treatment group in annotation_colors
annotation_colors <- list(
  Treatment = c("High" = "#E69F00", "Low" = "#009E73")
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
         main = "Z-Score Normalized logCPM for genes with Treatment Grouping")





# Create treatment labels for annotation, matching custom_column_labels order
treatment_labels <- c("High", "High", "High", "High", "High", "High", "High", "High",
                      "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", 
                      "High", "High", "High", "High", "High", "High", "High", "High",
                      "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", 
                      "Low", "Low", "Low", "Low", "Low", "Low", "Low", "High", "High",
                      "High", "High", "Low", "Low", "Low", "Low", "High", "High", "High",
                      "High", "High", "High", "High", "High", "Low", "Low", "Low", "Low",
                      "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low",
                      "Low", "Low", "High", "High", "High", "High", "High", "High",
                      "High", "High")

# Define the annotation data frame with row names matching custom_column_labels
annotation_col <- data.frame(Treatment = factor(treatment_labels, levels = c("High", "Low")))
rownames(annotation_col) <- custom_column_labels

# Confirm that `annotation_col` is properly structured
print(dim(annotation_col))  # Should show 80 rows and 1 column
print(head(annotation_col))  # Preview to ensure it looks correct

# Define colors for each treatment group in annotation_colors
annotation_colors <- list(
  Treatment = c("High" = "#E69F00", "Low" = "#009E73")
)

# Reorder the columns in `heatmap_data_filtered` and `annotation_col` based on treatment labels
ordered_indices <- order(annotation_col$Treatment)
heatmap_data_ordered <- heatmap_data_filtered[, ordered_indices]
annotation_col_ordered <- annotation_col[ordered_indices, , drop = FALSE]  # Ensure it remains a data frame

# Reorder the custom column labels to match the new order
custom_column_labels_ordered <- custom_column_labels[ordered_indices]

# Debugging Step: Check if the dimensions align after ordering
print(dim(heatmap_data_ordered))            # Print dimensions of the heatmap data
print(dim(annotation_col_ordered))          # Print dimensions of the reordered annotation
print(identical(colnames(heatmap_data_ordered), rownames(annotation_col_ordered))) # Should be TRUE

# Plot the heatmap with reordered data and annotations, without clustering columns
pheatmap(heatmap_data_ordered, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # Do not cluster columns; keep original order
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontface = "bold",
         scale = "none",  # Data is already z-scored
         color = inferno(256),
         labels_col = custom_column_labels_ordered, 
         annotation_col = annotation_col_ordered,  # Add the reordered column annotations for treatments
         annotation_colors = annotation_colors,  # Apply colors for each treatment group
         annotation_legend = TRUE,  # Show legend for the annotations
         main = "Z-Score Normalized logCPM dre-miR-10b-5p",
         # Adjust font size and row height
         fontsize_row = 3.5,         # Set font size for row labels
         fontsize_col = 6,         # Set font size for column labels
         cellheight = 3,          # Set row height
         cellwidth = 5)           # Optional: Set column width if you need to resize columns as well
