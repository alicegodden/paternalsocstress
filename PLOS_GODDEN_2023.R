# Second R scripts- PLOS GODDEN 2023.R

library(DESeq2)
library(tidyverse)
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
BiocManager::install('pheatmap')
library(ggplot2)

#use for rRNA seq and sub in the right file names 
countData <- read.csv('Supp_File_11_RawReadCounts_PIRNAS_6_matchedPairs.xlsx')
head(countData)


#delete first column of row numbers

df = subset(countData, select = -c( 1) )
rownames(df) <- countData[ , 1]
head(df)

mode(df)

colData <- read.csv('Supp_File_10_40_mirna_metadata_LF.csv')
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
                              design = ~ Treatment) 
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
dds$condition <-relevel(dds$condition, ref = "Control") #control
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
write.csv(normalized_counts, file="normalised_counts_mirnas.csv")
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
plotPCA(ldat, intgroup="Treatment") +
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
topVarGenes <- names(sort(row_variances, decreasing = TRUE)[1:30])

# Extract data for the top variable genes
heatmap_data <- norm_counts[topVarGenes, ]


# Plot the heatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontface = "bold",
         scale = "row",  # Normalize rows (genes) to have mean = 0 and sd = 1
         color= inferno(256),
         annotation_col = as.data.frame(colData(dds)[, "Treatment", drop=FALSE]))


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
