# Title: R scripts for paternal social stress paper
# Author : Dr. Alice M. Godden
library(DESeq2)
library(tidyverse)
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
BiocManager::install('pheatmap')
library(ggplot2)

#use for rRNA seq and sub in the right file names 
countData <- read.csv("mir_matched_counts_200.csv")
head(countData)


#delete first column of row numbers

df = subset(countData, select = -c( 1) )
rownames(df) <- countData[ , 1]
head(df)

mode(df)

colData <- read.csv('willian_metadata_ordered.csv')
head(colData)
#make first column of colData into row names
rownames(colData) <- colData [,1]  
head(colData)
colData = subset(colData, select = -c(1) )
head(colData)

#coldata$Treatment <- as.factor(coldata$Treatment)
#coldata$Reads <- as.factor(coldata$Reads)

#making sure the row names in colData match column names in df (your countdata)
all(colnames(df) %in% rownames(colData))
# if you get true then this is correct

#to ensure they are in same order
all(colnames(df) == rownames(colData))
# if true then this is correct


#one factor- miRNA and piRNA
dds <- DESeqDataSetFromMatrix(countData=round (df),
                              colData = colData,
                              design = ~ Group)

# for piRNA only DESeq to estimate size factors
#dds <- estimateSizeFactors(dds)
# Estimate gene-wise dispersions
#dds <- estimateDispersionsGeneEst(dds)
# Set the gene-wise dispersion estimates as the final dispersion estimates
#dispersions(dds) <- mcols(dds)$dispGeneEst
# Continue with the DESeq workflow
#dds <- nbinomWaldTest(dds)


#for two factor mRNA analysis to include female id as a factor
dds <- DESeqDataSetFromMatrix(countData=round (df),
                              colData = colData,
                              design = ~ condition+Female_id) 


#Run DEseq.
dds<-DESeq(dds)

dds


#pre-filtering: removing rows with low gene counts
#keeping rows that have at least 10 reads total


keep <- rowSums(counts(dds)) >= 10

dds <- dds [keep,]

dds

 
#step3 run DESeq
dds <- DESeq (dds)
resultsNames(dds)
dds$condition <-relevel(dds$condition, ref = "Low") #control
res <- results(dds)

#results(dds, contrast=list(c("Treatment_Temperature_vs_Control", "TreatmentTemperature.ReadsGood"  )))

res

#Explore results
summary(res)


res0.05 <- results(dds, alpha =0.05)
summary(res0.05)

write.csv(res, file="file.csv")

#to get normalised counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file="normalised_counts_mir10b.csv")
# contrasts

resultsNames(dds)



#PCA EXPLORER used to make bar plots for genes in mRNA figure
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pcaExplorer")

browseVignettes("pcaExplorer")
library("pcaExplorer")
library(ggplot2)

pcaExplorer(dds = dds)

#plotting own pca plots

nsub=nrow(dds)
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup = "Treatment", ntop = 300, returnData = TRUE)
pcaData <- plotPCA(rld, intgroup=c("Treatment"), returnData=TRUE)



###Final pca scripts
pcaData <- plotPCA(rld, intgroup=c("condition", "Female_id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=Female_id, shape=condition)) +
  geom_point(size=6) +
  scale_color_manual(values = c("Navy", "maroon4", "orchid", "plum", "Blue", "Pink", "Violet", "Purple", "Grey", "Black")) +
  scale_shape_manual(values=seq(0, 15)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +  # Set the background to white
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=16),
        
        #panel.grid.major = element_blank(),  # Remove grid lines
        #panel.grid.minor = element_blank(),
        axis.line = element_line(),  # Add axis lines
        axis.ticks = element_line())  # Add tick marks
      

-----

#pca used for Figure 1
ggplot(pcaData, aes(PC1, PC2, fill = Female_id, shape = condition, color = Female_id)) +
  geom_point(size = 6, alpha = 1) +
  scale_fill_manual(values = c("Navy", "maroon4", "orchid", "plum", "Blue", "Pink", "Violet", "Purple", "Grey", "Black")) +
  scale_shape_manual(values = c(21, 22, 21, 22, 21, 22, 21, 22, 21, 22)) +
  scale_color_manual(values = c("Navy", "maroon4", "orchid", "plum", "Blue", "Pink", "Violet", "Purple", "Grey", "Black")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.line = element_line(),
    axis.ticks = element_line()
  )


  #option 2 didn't use this but makes pca points hollow- no colour fill just outline
ggplot(pcaData, aes(PC1, PC2, color = Female_id, shape = condition)) +
  geom_point(size = 6) +
  scale_color_manual(values = c("Navy", "maroon4", "orchid", "plum", "Blue", "Pink", "Violet", "Purple", "Grey", "Black")) +
  scale_shape_manual(values = seq(0, 15)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +  # Set the background to white
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.line = element_line(),  # Add axis lines
    axis.ticks = element_line()  # Add tick marks
  )

########## Making mirna and piRNA pca plot
rld <- rlog(dds, blind=FALSE)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
nbinomLRT(dds, reduced = ~1)


####################### to then make your pca plot for piRNA and miRNA
ldat <- normTransform(dds)

ldat <- normTransform(dds, f = log2, pc =1)

rld <-rlog(ldat)
plotPCA(ldat, intgroup="Group") +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title = element_text(size = 15))

counts(ldat)



rld <-rlog(ldat)
plotPCA(ldat, intgroup="Treatment") +
  geom_point(size=6) +
  scale_color_manual(values = c("Pink", "Blue", "Purple", "Violet")) +
  scale_shape_manual(values=seq(0, 15)) +
  theme_bw() +  # Set the background to white
  theme(axis.text.x = element_text(size = 16, face='bold', colour = 'black'),
        axis.text.y = element_text(size = 16, face='bold', colour = 'black'),
        axis.title = element_text(size = 18, face='bold', colour = 'black'),
        legend.title = element_text(size = 18, face='bold', colour = 'black'),
        legend.text = element_text(size=16, face='bold', colour = 'black'),
        
        #panel.grid.major = element_blank(),  # Remove grid lines
        #panel.grid.minor = element_blank(),
        axis.line = element_line(),  # Add axis lines
        axis.ticks = element_line())  # Add tick marks

##################

#Other pca plots, same scheme but are one factor and less colours

ggplot(pcaData, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=6) +
  scale_color_manual(values = c("Pink", "Blue", "Purple", "Violet")) +
  scale_shape_manual(values=seq(0, 15)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +  # Set the background to white
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size=16),
        
        #panel.grid.major = element_blank(),  # Remove grid lines
        #panel.grid.minor = element_blank(),
        axis.line = element_line(),  # Add axis lines
        axis.ticks = element_line())  # Add tick marks

###Plotting normalised counts heatmaps for sRNAs
# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(viridis)
# Normalize and log-transform the data with a pseudo count of 1
ldat <- normTransform(dds, f = log2, pc = 1)

# Extract normalized counts
norm_counts <- assay(ldat)

# Calculate row variances, explicitly setting useNames = TRUE
row_variances <- rowVars(norm_counts, useNames = TRUE)

# Order the genes by variance and select the top 30 most variable genes
topVarGenes <- names(sort(row_variances, decreasing = TRUE)[1:46])

# Extract data for the top variable genes
heatmap_data <- norm_counts[topVarGenes, ]


# Plot the heatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         main = "dre-miR-200a-5p",
         fontsize_row = 4,
         fontsize_col = 8,
         cellheight = 4,
         cellwidth = 8,
         fontface = "bold",
         scale = "none",  # Normalize rows (genes) to have mean = 0 and sd = 1
         color= inferno(256),
         annotation_col = as.data.frame(colData(dds)[, "condition", drop=FALSE]))

# grouping High and Low conditions together
library(pheatmap)


# Example data (replace with your actual heatmap_data and dds)
heatmap_data <- matrix(rnorm(100), nrow = 10)  # Example heatmap data
dds <- data.frame(condition = rep(c("Low", "High"), each = 5))  # Example condition data

# Specify the desired column order
desired_order <- c(
  "E100.1", "E100.2", "E100.3", "E100.4",
  "E101.1", "E101.2", "E101.3", "E101.4",
  "E104.1", "E104.2", "E104.3", "E104.4",
  "E105.1", "E105.2", "E105.3", "E105.4",
  "E110.1", "E110.2", "E110.3", "E110.4",
  "E112.4", "E112.5", "E112.7", "E112.8",
  "E113.10", "E113.11", "E113.13", "E113.9",
  "E82.1", "E82.2", "E82.3", "E82.5",
  "E83.1", "E83.3", "E83.4", "E83.5",
  "E102.1", "E102.2", "E102.3", "E102.4",
  "E103.1", "E103.2", "E103.3", "E103.4",
  "E106.1", "E106.2", "E106.3", "E106.4",
  "E107.1", "E107.2", "E107.3", "E107.4",
  "E108.1", "E108.2", "E108.3", "E108.4",
  "E109.1", "E109.2", "E109.3", "E109.4",
  "E111.1", "E111.2", "E111.3", "E111.4",
  "E114.1", "E114.3", "E114.4", "E114.5",
  "E115.1", "E115.2", "E115.4", "E115.5",
  "E80.1", "E80.2", "E80.3", "E80.4",
  "E81.1", "E81.2", "E81.4", "E81.5"
)

# Ensure heatmap_data and dds are ordered according to desired_order
order <- match(colnames(heatmap_data), desired_order)
heatmap_data_ordered <- heatmap_data[, order]

# Create the heatmap
pheatmap(heatmap_data_ordered,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Heatmap with Specific Column Order",
         font.face = "bold",
         scale = "none",  # Normalize rows (genes) to have mean = 0 and sd = 1
         color = inferno(256),
         annotation_col = as.data.frame(dds),
         annotation_colors = list(condition = c("Low" = "blue", "High" = "red")))

## To specify differentially expressed sRNAs-
# Sample heatmap data (replace this with your actual data)
# Define the genes you want to display

# Load necessary libraries
library(DESeq2)   # Assuming norm_counts is derived from DESeq2
library(pheatmap)
library(viridis)  # For the inferno color palette

# Define the normalized counts data
# Assuming norm_counts is your normalized counts matrix
# norm_counts <- assay(dds)  # Uncomment and set this if needed

# Define the genes you want to display
genes_to_display <- c("dre-miR-129-5p", "dre-miR-184", "dre-miR-181a-5-3p", "dre-miR-183-5p",
                      "dre-miR-193b-5p", "dre-miR-10b-5p", "dre-miR-200a-5p")

# Check which genes are available in the dataset
available_genes <- rownames(norm_counts)
selected_genes <- genes_to_display[genes_to_display %in% available_genes]

# Subset the norm_counts to include only the specified genes and all samples
subset_norm_counts <- norm_counts[selected_genes, , drop = FALSE]

# Check if there are at least 2 rows after subsetting
if (nrow(subset_norm_counts) < 2) {
  stop("Insufficient rows (genes) to perform clustering. The following genes are available: ", paste(selected_genes, collapse = ", "))
}

# Ensure annotation_col matches the columns of subset_norm_counts
annotation_col <- as.data.frame(colData(dds)[, "Treatment", drop=FALSE])
annotation_col <- annotation_col[colnames(subset_norm_counts), , drop = FALSE]


# Define custom labels for columns, give them meaningful names
custom_column_labels <- c("High.2_1", "High.2_2", "High.2_3", "High.2_4", "High.2_5", "High.2_6", "High.2_7", "High.2_8", "High.2_9", "High.2_10",
                          "Low.1_1", "Low.1_2", "Low.1_3", "Low.1_4", "Low.1_5", "Low.1_6", "Low.1_7", "Low.1_8", "Low.1_9", "Low.1_10",
                          "Low.2_1", "Low.2_2", "Low.2_3", "Low.2_4", "Low.2_5", "Low.2_6", "Low.2_7", "Low.2_8", "Low.2_9", "Low.2_10",
                          "High.1_1", "High.1_2", "High.1_3", "High.1_4", "High.1_5", "High.1_6", "High.1_7", "High.1_8", "High.1_9", "High.1_10")

# Check that the length of custom_column_labels matches the number of columns in subset_norm_counts
if (length(custom_column_labels) != ncol(subset_norm_counts)) {
  stop("Length of custom_column_labels does not match the number of columns in subset_norm_counts.")
}

# Generate the heatmap
pheatmap(subset_norm_counts,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "none",  # No scaling
         color = inferno(256),
         labels_col = custom_column_labels, 
         annotation_col = annotation_col,
         fontsize_row = 12,  # Set font size for row names
         fontsize_col = 12,  # Set font size for column names
         cellwidth = 9.5,
         cellheight = 20,
         fontface = "bold")  # Bold text

# to plot those with > 100 normalised reads 

# Normalize and log-transform the data with a pseudo count of 1
ldat <- normTransform(dds, f = log2, pc = 1)

# Extract normalized counts
norm_counts <- assay(ldat)

# Filter rows where all counts are greater than 100
filtered_counts <- norm_counts[apply(norm_counts, 1, function(x) all(x > 20)), ]

# Calculate row variances, explicitly setting useNames = TRUE
row_variances <- rowVars(filtered_counts, useNames = TRUE)

# Order the genes by variance and select the top 30 most variable genes
topVarGenes <- names(sort(row_variances, decreasing = TRUE)[1:30])

# Extract data for the top variable genes
heatmap_data <- filtered_counts[topVarGenes, ]

# Plot the heatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         labels_col = custom_column_labels,
         fontface = "bold",
         scale = "none",  # Do not scale rows (genes)
         color = inferno(256),
         annotation_col = as.data.frame(colData(dds)[, "Treatment", drop = FALSE]))


# for pirna as they have long row names

# Normalize and log-transform the data with a pseudo count of 1
ldat <- normTransform(dds, f = log2, pc = 1)

# Extract normalized counts
norm_counts <- assay(ldat)

# Calculate row variances, explicitly setting useNames = TRUE
row_variances <- rowVars(norm_counts, useNames = TRUE)

# Order the genes by variance and select the top 30 most variable genes
topVarGenes <- names(sort(row_variances, decreasing = TRUE)[1:30])

# Extract data for the top variable genes
heatmap_data <- norm_counts[topVarGenes, ]

# Crop row names to "NC_" or "CHROM" plus the next 30 characters
rownames(heatmap_data) <- sub("{0,30}", "\\1", rownames(heatmap_data))

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


#now we want the genes in our deseq2 list to be matched to ensembl names to gene names
library(biomaRt)
library(tidyverse)
library(org.Dr.eg.db)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Dr.eg.db")

#read in your significant gene ids
ensembl.ids <- read.delim('genes_id_AMG_te.csv', header =F)
head (ensembl.ids)

#biomart
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
datasets

ensembl.con <- useMart("ensembl", dataset = "drerio_gene_ensembl")

attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

IDs <- getBM(attributes = c('ensembl_gene_id', 'gene_biotype', 'external_gene_name'),
      filters = "ensembl_gene_id",
      values = ensembl.ids$V1,
      mart = ensembl.con)


head(IDs)
write.csv(IDs, file="genes_AMG_te.csv")



#load in your deseq2 output file, and add ensembl_id to first column header
resid <- read.csv('te_tocontrol_deseq2.csv')
# Merging two tables, syntax showed here are in full forms
annot.table <- merge(x = resid, y = IDs,
                     by.x = "ensembl_id", by.y =
                       "ensembl_gene_id", all.x = T, all.y = F )
# export your annotated results
write.csv(annot.table, file="annotated_te_tocontrol_deseq2.csv")


#now make a heatmap


#Enhanced volcano
library (tidyverse)
library(EnhancedVolcano)

# assign the results file name to a variable
deseq_results_file <- 'volcano.tsv'

# load data
deseq_results <-read_tsv(deseq_results_file,
           col_types = cols(Chr = col_character(),
                            Strand = col_character()))


head(deseq_results)
glimpse(deseq_results)
View(deseq_results)

EnhancedVolcano(
  deseq_results, # results data frame
  lab = deseq_results$TE,
  x = 'log2FoldChange', # column name of log2 fold change
  y = 'padj' # column name of adjusted pvalue
)


#to add the top 10 most sig de genes by p calue
genes_to_label <- deseq_results %>%
  arrange(padj) %>%#change log2fc to padj for by sig p value
  pull(miRNA) %>%
  head (30)

( ?EnhancedVolcano )

EnhancedVolcano(
  deseq_results, # results data frame
  lab = deseq_results$miRNA,
  x = 'log2FoldChange', # column name of log2 fold change
  y = 'padj', # column name of adjusted pvalue
  col = c("grey", "darkgreen", "navy", "gold"), 
  pCutoff = 50e-03,
  FCcutoff = log2(2.0),
  ylim = c(0,2.5),
  xlim = c(-3,3),
  min.segment.length = 0.1,
  labFace = "bold",
  labSize = 2.5,
  pointSize = 1.5,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  title = 'Enhanced Volcano',
  subtitle = 'sRNA-seq: Testes hairpin miRNAs most differentially expressed by padj'
)


library(tidyverse)
#Heatmap
# assign the results file name to a variable
deseq_results_file <- ('VOLCANO_na_annotated_ov_tocontrol_deseq2.tsv')
# load data
deseq_results <-read_tsv(deseq_results_file,
                         col_types = cols(Chr = col_character(),
                                          Strand = col_character()))



#heatmap
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
library(readr)
#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

data <- read.csv("piRNA_pheatmap.csv")#, comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names

head(mat_data)
#########################################################
### C) Customizing and plotting the heat map
#########################################################
#Spare heatmap script
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green

usethis::edit_r_environ("project") 

#spare heatmap script
heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering
library(pheatmap)



browseVignettes("DESeq")
pheatmap(countData)
sig_counts <- deseq_results 
countData = as.matrix(countData)

#make first column of colData into row names
rownames(countData) <- countData [,1]  
head(countData)
countData = subset(countData, select = -c(1) )
head(countData)

#USed the pheatmap script

pheatmap(mat_data)

pheatmap(mat_data,
         scale = "row")

library(viridis)
scales::show_col(inferno(10000), cex_label = 0.6)

scales::show_col(colorRampPalette(inferno(10))(10000),
                 cex_label = 0.6)
pheatmap(mat_data,
         scale = "row", 
         show_rownames = TRUE,
         angle_col = 45,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         treeheight_row = 100,
         fontsize_row = 6,
         color = colorRampPalette(inferno(10))(10000))
