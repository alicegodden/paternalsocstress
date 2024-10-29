# Hypergeometric test of miRNA target genes compared to all DE genes
# Parameters
N <- 15198  # total number of genes in your background
M <- 358    # total number of miRNA target genes in the background
n <- 612    # number of DE genes
k <- 50     # number of DE genes that are also miRNA targets (replace this with the actual overlap count)

# Hypergeometric test
p_value <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)

# Print the p-value to interpret enrichment
print(p_value)
