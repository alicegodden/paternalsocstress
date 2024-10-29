
# 15198 is total number of genes 
# 612 is number of significantly differentially expressed genes
# 515 is total number of genes targeted by our mirnas
# 358 is number of significanly differentially expressed genes targeted by our miRNAs
#The table deTable now has the following format, matching the structure in the example:
  
#  In Gene Set (miRNA Target)	Out of Gene Set (Non-miRNA Target)
#DE (yes)	358	254
#Not DE (no)	157	14429
# Construct the contingency table in the same format as the example
deTable <- matrix(c(358, 254, 157, 14429),
                  nrow = 2,
                  dimnames = list(DE = c("yes", "no"),
                                  GeneSet = c("in", "out")))

# Display the table
deTable

# Perform Fisher's exact test
fisher.test(deTable, alternative = "greater")


# Parameters
N <- 15198      # Total genes
M <- 515        # Total miRNA target genes in the gene set
n <- 612        # Total DE genes
k <- 358        # DE genes that are also in the gene set (miRNA targets)

# Hypergeometric test
p_value <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)

# Print the p-value
print(paste("P-value:", p_value))

