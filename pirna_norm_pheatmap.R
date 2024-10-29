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
  "NC_007115.7:43701007-43707025",
  "12225229-12234949",
  "NC_007133.7:32643430-32648896",
  "NC_007115.7:52720448-52740603",
  "32505005-32521009",
  "24796020-24822028",
  "NC_007123.7:14899002-14913024",
  "747001-760987",
  "57940016-57946585",
  "14845000-14858022",
  "63048-99933",
  "NC_007118.7:2885044-2893968",
  "7007002-37018889",
  "15543003-15555996",
  "29825000-29835935",
  "29842000-29853219",
  "12939021-12949958",
  "12957011-12967948",
  "60877023-60885024",
  "8206091-8220020",
  "NC_007125.7:13384000-13389778",
  "13797593-13804958",
  "27056031-27069027",
  "13181000-13191573",
  "16363000-16370764",
  "NC_007115.7:41646063-41653819",
  "NC_007131.7:15827011-15836010",
  "NC_007128.7:40817625-40825762",
  "57283001-57302029",
  "NC_007129.7:33974002-33981018",
  "57344003-57355026",
  "2325020-2330026",
  "NC_007115.7:73175084-73181996",
  "52547012-52558402",
  "15130054-15141020",
  "56143038-56151691",
  "56159035-56167792",
  "17874034-17886001",
  "NC_007116.7:47663015-47686016",
  "15825004-15836000",
  "35696144-35707638",
  "NC_007115.7:46831029-46839861",
  "43507012-43526027",
  "5785001-5802968",
  "8068108-8076992",
  "3976000-23991019",
  "3654000-13670019",
  "32121033-32144977",
  "13021007-13051030",
  "NC_007129.7:32561068-32568584",
  "47728027-47753986",
  "23686113-23695941",
  "11251294-11261493",
  "73380021-73389987",
  "NC_007131.7:17034025-17042856",
  "10080045-10093028")

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
  "4_DNA:CHROMOSOME_CHROMOSOME:GRCZ11:4:1:78093715:1_REF_27056031-27069027_(+-0_BP)_DIRECTIONALITY:_MONO:PLUS",
  "CHR_ALT_CTG4_1_12_DNA:CHROMOSOME_CHROMOSOME:GRCZ11:CHR_ALT_CTG4_1_12:1:78082549:1_HAP_27056031-27069027_(+-0_BP)_DIRECTIONALITY:_MONO:PLUS",
  "CHR_ALT_CTG18_1_17_DNA:CHROMOSOME_CHROMOSOME:GRCZ11:CHR_ALT_CTG18_1_17:1:51050792:1_HAP_16363000-16370764_(+-0_BP)_DIRECTIONALITY:_MONO:PLUS"
)
heatmap_data_filtered <- heatmap_data[!(rownames(heatmap_data) %in% rows_to_remove), ]

# Define custom column labels for the heatmap
custom_column_labels <- c("High.2_2", "High.2_6", "High.2_7", "High.2_8", "High.2_9", "High.2_10",
                          "Low.1_4", "Low.1_5", "Low.1_6", "Low.1_7", "Low.1_9", "Low.1_10",
                          "Low.2_2", "Low.2_3", "Low.2_4", "Low.2_5", "Low.2_6", "Low.2_9",
                          "High.1_1", "High.1_3", "High.1_4", "High.1_6", "High.1_7", "High.1_10")

# Verify that the column names in heatmap_data_filtered match custom_column_labels
colnames(heatmap_data_filtered) <- custom_column_labels  # Ensure columns have custom names

# Create treatment labels for annotation, matching custom_column_labels order
treatment_labels <- c("High2", "High2", "High2", "High2", "High2", "High2",
                      "Low1", "Low1", "Low1", "Low1", "Low1", "Low1",
                      "Low2", "Low2", "Low2", "Low2", "Low2", "Low2",
                      "High1", "High1", "High1", "High1", "High1", "High1")

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
         main = "Heatmap of Z-Score Normalized logCPM for Selected piRNAs with Treatment Grouping")

# Create treatment labels for annotation
treatment_labels <- c("High2", "High2", "High2", "High2", "High2", "High2",
                      "Low1", "Low1", "Low1", "Low1", "Low1", "Low1",
                      "Low2", "Low2", "Low2", "Low2", "Low2", "Low2",
                      "High1", "High1", "High1", "High1", "High1", "High1")

# Annotation data for columns
annotation_col <- data.frame(Treatment = factor(treatment_labels))
rownames(annotation_col) <- custom_column_labels

# Define colors for each treatment group
annotation_colors <- list(
  Treatment = c("High1" = "#E69F00", "High2" = "#56B4E9", "Low1" = "#009E73", "Low2" = "#F0E442")
)

## Plot the heatmap
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
         main = "Z-Score Normalized logCPM for significantly differentially expressed piRNA clusters")
