# Load necessary libraries
library(stringr)
library(readr)

# Make sure you already have your dds object from the main R script

# Replace with your actual file path
counts_file <- "pirna_24_rawcounts_use.csv"

# Read the counts data from CSV file
counts_data <- read.csv(counts_file, row.names = 1)

# List of row names
# Normalize and log-transform the data with a pseudo count of 1
ldat <- normTransform(dds, f = log2, pc = 1)

# Extract normalized counts
norm_counts <- assay(ldat)

# Calculate row variances, explicitly setting useNames = TRUE
row_variances <- rowVars(norm_counts, useNames = TRUE)

# List of row names to match
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
  matched_row <- rownames(norm_counts)[grep(row_name, rownames(norm_counts))]
  return(matched_row)
}

# Match each row name in the list
matched_rows <- sapply(row_names_to_match, match_row)


# Extract data for the matched rows
heatmap_data <- norm_counts[unlist(matched_rows), ]


# Show matching rows in norm_counts
matching_norm_counts <- norm_counts[heatmap_rows, ]
# Set rownames of heatmap_data to row_names_to_match
rownames(heatmap_data) <- row_names_to_match

# Modify row names to add "CHROM:" if row name is numeric
new_row_names <- row_names_to_match
numeric_pattern <- "^[0-9-]+$"
for (i in seq_along(new_row_names)) {
  if (grepl(numeric_pattern, new_row_names[i])) {
    new_row_names[i] <- paste("CHROM:", new_row_names[i], sep = "")
  }
}

# Set rownames of heatmap_data
rownames(heatmap_data) <- new_row_names


# Plot the heatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontface = "bold",
         scale = "none",  # Normalize rows (genes) to have mean = 0 and sd = 1
         color = inferno(256),
         annotation_col = as.data.frame(colData(dds)[, "Group", drop = FALSE]))

# to give column names meaningful names to match other plots
custom_column_labels <- c("High.2_2", "High.2_6", "High.2_7", "High.2_8", "High.2_9", "High.2_10",
                          "Low.1_4", "Low.1_5", "Low.1_6", "Low.1_7", "Low.1_9", "Low.1_10",
                          "Low.2_2", "Low.2_3", "Low.2_4", "Low.2_5", "Low.2_6", "Low.2_9",
                          "High.1_1", "High.1_3", "High.1_4", "High.1_6", "High.1_7", "High.1_10")


# Plot the heatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontface = "bold",
         scale = "none",  # Normalize rows (genes) to have mean = 0 and sd = 1
         color = inferno(256),
         labels_col = custom_column_labels, 
         annotation_col = as.data.frame(colData(dds)[, "Group", drop = FALSE]))

# to identify and remove duplicate rows
# Assuming norm_counts is already defined and contains your data
# Define row_names_to_match as before

# Function to match row names
match_row <- function(row_name) {
  matched_rows <- rownames(norm_counts)[grep(row_name, rownames(norm_counts))]
  return(matched_rows)
}

# Match each row name in the list
matched_rows <- sapply(row_names_to_match, match_row, simplify = FALSE)

# Identify rows that match to more than one row in norm_counts
multiple_matches <- matched_rows[sapply(matched_rows, length) > 1]

# Extract data for the matched rows
heatmap_data_multiple <- norm_counts[unlist(multiple_matches), ]

# Print row names that have multiple matches
cat("Row names with multiple matches:\n")
print(names(multiple_matches))

# Print the matched rows in norm_counts
cat("\nMatching rows in norm_counts:\n")
print(heatmap_data_multiple)


# checked against full names need to remove:
#4_DNA:CHROMOSOME_CHROMOSOME:GRCZ11:4:1:78093715:1_REF_27056031-27069027_(+-0_BP)_DIRECTIONALITY:_MONO:PLUS
#CHR_ALT_CTG4_1_12_DNA:CHROMOSOME_CHROMOSOME:GRCZ11:CHR_ALT_CTG4_1_12:1:78082549:1_HAP_27056031-27069027_(+-0_BP)_DIRECTIONALITY:_MONO:PLUS 
#CHR_ALT_CTG18_1_17_DNA:CHROMOSOME_CHROMOSOME:GRCZ11:CHR_ALT_CTG18_1_17:1:51050792:1_HAP_16363000-16370764_(+-0_BP)_DIRECTIONALITY:_MONO:PLUS     
#from heatmap_data

# Rows to remove
rows_to_remove <- c(
  "4_DNA:CHROMOSOME_CHROMOSOME:GRCZ11:4:1:78093715:1_REF_27056031-27069027_(+-0_BP)_DIRECTIONALITY:_MONO:PLUS",
  "CHR_ALT_CTG4_1_12_DNA:CHROMOSOME_CHROMOSOME:GRCZ11:CHR_ALT_CTG4_1_12:1:78082549:1_HAP_27056031-27069027_(+-0_BP)_DIRECTIONALITY:_MONO:PLUS",
  "CHR_ALT_CTG18_1_17_DNA:CHROMOSOME_CHROMOSOME:GRCZ11:CHR_ALT_CTG18_1_17:1:51050792:1_HAP_16363000-16370764_(+-0_BP)_DIRECTIONALITY:_MONO:PLUS"
)

# Filter out the rows from heatmap_data
heatmap_data_filtered <- heatmap_data[!(rownames(heatmap_data) %in% rows_to_remove), ]

# Plot the heatmap with filtered data
pheatmap(heatmap_data_filtered, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontface = "bold",
         scale = "none",  # Normalize rows (genes) to have mean = 0 and sd = 1
         color = inferno(256),
         annotation_col = as.data.frame(colData(dds)[, "Group", drop = FALSE]))

#go back use above script after running this:
heatmap_data <- heatmap_data_filtered



