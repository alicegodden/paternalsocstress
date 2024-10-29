# Read in the list of DE genes
de_genes_df <- read.csv("DE_GENES.csv", header = FALSE)
DE_genes <- de_genes_df[[1]]  # Assuming gene IDs are in the first column

# Read in the list of miRNA target genes
mirna_targets_df <- read.csv("genes_targeted_mirna.csv", header = FALSE)
miRNA_targets <- mirna_targets_df[[1]]  # Assuming gene IDs are in the first column


k <- length(intersect(DE_genes, miRNA_targets))


# Hypergeometric test of miRNA target genes compared to all DE genes
# Parameters
N <- 15198  # total number of genes in your background
M <- 358    # total number of miRNA target genes in the background
n <- 612    # number of DE genes


# Hypergeometric test
p_value <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)

# Print the p-value to interpret enrichment
print(p_value)


# Calculate Expected Overlap
expected_overlap <- (M / N) * n

# Calculate Fold Enrichment
fold_enrichment <- k / expected_overlap

# Print the results
print(paste("P-value:", p_value))
print(paste("Fold Enrichment:", fold_enrichment))
