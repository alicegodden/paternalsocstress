# Title: R scripts used in Paternal social stress data analysis
# Author: Dr. Alice M. Godden

library(DESeq2)
library(tidyverse)
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
BiocManager::install('pheatmap')

#use for rRNA seq and sub in the right file names 
countData <- read.csv('Processed_Data_40_idep.csv')
head(countData)


#delete first column of row numbers

df = subset(countData, select = -c( 1) )
rownames(df) <- countData[ , 1]
head(df)

#make first column of colData into row names
#rownames(countData) <- countData [,1]  
#head(countData)
#countData = subset(countData, select = -c(1) )
#head(countData)

mode(df)



colData <- read.csv('40_mirna_metadata_LF.csv')
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


#one factor- to avoid this error: Error in DESeqDataSet(se, design = design, ignoreRank) : 
#some values in assay are not integers
#dds <- DESeqDataSetFromMatrix(countData=round (df),

dds <- DESeqDataSetFromMatrix(countData=round (df),
                              colData = colData,
                              design = ~ Treatment) 

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


#to get normalised counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file="willian_rna_norm_counts_deseq2.csv")

res0.05 <- results(dds, alpha =0.05)
summary(res0.05)

write.csv(res, file="female_telescope_deseq2.csv")

#PCA plotting
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup = "condition", "Female_id", ntop = 300, returnData = TRUE)
pcaData <- plotPCA(rld, intgroup=c("condition", "Female_id"), returnData=TRUE)


pcaData <- plotPCA(rld, intgroup=c("condition", "Female_id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=group, shape=age)) +
  geom_point(size=6) +
  scale_color_manual(values = c("Navy", "maroon4", "orchid", "plum", "Blue", "Pink", "Violet", "Purple", "Grey", "Black")) +
  scale_shape_manual(values=seq(0, 15)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +  # Set the background to white
  theme(axis.text.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text.y = element_text(size = 16, color = "black", face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size=16, face = "bold"),
        
        #panel.grid.major = element_blank(),  # Remove grid lines
        #panel.grid.minor = element_blank(),
        axis.line = element_line(),  # Add axis lines
        axis.ticks = element_line())  # Add tick marks

#pca for miRNA and piRNA

########## Making piRNA pca plot
rld <- rlog(dds, blind=FALSE)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
nbinomLRT(dds, reduced = ~1)


####################### to then make your pca plot for piRNA
ldat <- normTransform(dds)

ldat <- normTransform(dds, f = log2, pc =1)

rld <-rlog(ldat)
plotPCA(ldat, intgroup="Group") +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title = element_text(size = 15))

counts(ldat)



rld <-rlog(ldat)
plotPCA(rld, intgroup="Treatment") +
  geom_point(size=6) +
  scale_color_manual(values = c("Pink", "Blue", "Purple", "Violet")) +
  scale_shape_manual(values=seq(0, 15)) +
  theme_bw() +  # Set the background to white
  theme(axis.text.x = element_text(size = 16, face = 'bold', color= 'black'),
        axis.text.y = element_text(size = 16, face = 'bold', color = 'black'),
        axis.title = element_text(size = 18, face = 'bold'),
        legend.title = element_text(size = 18, face = 'bold'),
        legend.text = element_text(size=16, face = 'bold'),
        
        #panel.grid.major = element_blank(),  # Remove grid lines
        #panel.grid.minor = element_blank(),
        axis.line = element_line(),  # Add axis lines
        axis.ticks = element_line())  # Add tick marks

##################

#Other pca plots, same scheme but are one factor and less colours

ggplot(ldat, aes(PC1, PC2, color=Group)) +
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

#volcano plotting

# assign the results file name to a variable
deseq_results_file <- 'male_group_all_telescope_volcano.tsv'

# load data
deseq_results <-read_tsv(deseq_results_file,
                         col_types = cols(Chr = col_character(),
                                          Strand = col_character()))


head(deseq_results)
glimpse(deseq_results)
View(deseq_results)


library(EnhancedVolcano)
EnhancedVolcano(
  deseq_results, # results data frame
  lab = deseq_results$TE,
  x = 'log2FoldChange', # column name of log2 fold change
  y = 'padj' # column name of adjusted pvalue
)


#to add the top 10 most sig de genes by p calue
genes_to_label <- deseq_results %>%
  arrange(padj) %>%#change log2fc to padj for by sig p value
  pull(external_gene_name) %>%
  head (30)

( ?EnhancedVolcano )
# Install and load the Viridis package
install.packages("viridis")
library(viridis)

# Set the Viridis color palette
#set_palette("plasma")

# Load the EnhancedVolcano library
library(EnhancedVolcano)

# Your existing code for EnhancedVolcano
EnhancedVolcano(
  deseq_results, # results data frame
  lab = deseq_results$TE,
  x = 'log2FoldChange', # column name of log2 fold change
  y = 'padj', # column name of adjusted pvalue
  col = rocket(6), # Use Viridis color palette with 4 colors
  pCutoff = 50e-03,
  FCcutoff = log2(1.5),
  ylim = c(0, 8),
  xlim = c(-10, 10),
  min.segment.length = 0.1,
  labFace = "bold",
  labSize = 2.5,
  pointSize = 1.5,
  drawConnectors = TRUE,
  widthConnectors = 0.2,
  title = 'Enhanced Volcano',
  subtitle = 'Telescope Male All SEG v NEG TEs'
)


##heatmap


#normalised counts filter for bar chart
library(dplyr)
# Read the CSV file
data <- read.csv("tetrans_norm_counts_deseq2.csv", header = TRUE)
head(data)

# Find and print rows containing "hAT", "TDR", and "DNA" , not case sensitive will pick up lower and uppercase
filtered_data <- data[grep(("DNA|hAT|Tc1|TcMar|Harbinger|Enspm|Kolobok|Merlin|Crypton|PiggyBac|Dada|Zatar|Ginger|TDR|Polinton|Maverick|Acrobat|Looper|TZF|Angel|Mariner"), data$X, ignore.case = TRUE), ]
cat("Rows with 'DNA', 'hAT', 'Tc1', 'TcMar', 'Harbinger', 'Enspm', 'Kolobok', 'Merlin', 'Crypton', 'PiggyBac', 'Dada', 'Zatar', 'Ginger', 'TDR', 'Polinton', 'Maverick', 'Acrobat', 'Looper', 'TZF', 'Angel', or 'Mariner':\n")
print(filtered_data)

# Count and print the sum of each column (columns 2-13)
column_sums <- colSums(filtered_data[, 2:81])
cat("Sum of each column (columns 2-81):\n")
print(column_sums)
write.csv(column_sums, file="tetrans_norm_counts_dna.csv")



# bar chart for TEtranscripts
library(ggplot2)
library(tidyr)
library(ggpubr)

# Create a data frame from the image data
data <- data.frame(
  Treatment = c("High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low", "Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low"),
  gene = c("E100.1",	"E100.2",	"E100.3",	"E100.4",	"E101.1",	"E101.2",	"E101.3",	"E101.4",	"E104.1",	"E104.2",	"E104.3",	"E104.4",	"E105.1",	"E105.2",	"E105.3",	"E105.4",	"E110.1",	"E110.2",	"E110.3",	"E110.4",	"E112.4",	"E112.5",	"E112.7",	"E112.8",	"E113.1",	"E113.11",	"E113.13",	"E113.9",	"E82.1",	"E82.2",	"E82.3",	"E82.5",	"E83.1",	"E83.3",	"E83.4",	"E83.5",	"E102.1",	"E102.2",	"E102.3",	"E102.4",	"E103.1",	"E103.2",	"E103.3",	"E103.4",	"E106.1",	"E106.2",	"E106.3",	"E106.4",	"E107.1",	"E107.2",	"E107.3",	"E107.4",	"E108.1",	"E108.2",	"E108.3",	"E108.4",	"E109.1",	"E109.2",	"E109.3",	"E109.4",	"E111.1",	"E111.2",	"E111.3",	"E111.4",	"E114.1",	"E114.3",	"E114.4",	"E114.5",	"E115.1",	"E115.2",	"E115.4",	"E115.5",	"E80.1",	"E80.2",	"E80.3",	"E80.4",	"E81.1",	"E81.2",	"E81.4",	"E81.5"),
  LINE = c(172051.010343351,	181798.171233951,	167239.777767372,	192813.461754451,	234507.912857316,	198465.146952479,	219371.804003757,	254715.971981223,	177872.599878578,	195369.486516681,	186575.549064863,	201324.637555469,	179216.10968571,	186421.800763887,	185120.128684641,	198193.776832608,	193190.522875243,	193527.751991868,	181011.622541169,	184507.011935077,	226235.773212631,	197894.210778839,	189854.351853873,	184072.883442402,	174914.796878524,	172060.599557008,	171793.431070693,	168941.949434443,	180536.392262371,	169744.986519629,	174388.797247282,	186705.018684362,	177380.640450672,	185643.461091638,	173591.46622382,	182114.696361727,	231930.833119505,	219266.959940662,	205004.876206464,	183653.148517671,	202136.890583854,	207700.262207727,	195273.538685284,	189363.016369216,	196400.954065268,	186033.660168918,	184038.385310046,	186215.735106491,	184769.094221855,	183421.014720214,	184304.639842552,	180285.926030174,	194996.863373054,	185060.206594937,	232103.284254028,	213714.890062587,	201696.466066612,	215539.518364105,	192318.224264503,	191506.288687014,	208302.467984097,	196633.330097499,	189652.449112463,	198504.001771491,	214790.24181259,	187171.442263759,	179995.988899247,	189170.561691902,	175121.940174235,	167420.623928208,	180161.10058392,	197874.731506271,	184892.33691378,	177169.710953803,	181255.393670272,	184596.950712525,	182702.720566856,	178633.241812634,	182486.315808998,	174698.646682999),
  LTR = c(277780.208508856,	271662.579920513,	296765.936026165,	284761.0094366,	377902.905743328,	381430.290894118,	380749.370267615,	489972.54731661,	331035.299948178,	308977.340152634,	377095.568277416,	324245.090479094,	292806.353632439,	343807.728534855,	356311.348923484,	392302.309423938,	321100.436413208,	319298.147192926,	320190.350886054,	290906.780137328,	297408.829586074,	255594.055233358,	309829.121959588,	296220.638665861,	294594.997037001,	310879.856844734,	316555.937358353,	318786.886085994,	326885.371625091,	319104.167137172,	301124.746096668,	295943.277073925,	327896.431805257,	338734.92982271,	297646.111436376,	303651.843397811,	436543.626847255,	321077.594390848,	309852.342250699,	317763.019481653,	360952.993110389,	413393.838420643,	362947.175144834,	386651.88078083,	345211.088563561,	335475.845214327,	319404.718973342,	332128.615494627,	296204.415940727,	319459.834446284,	335160.67234884,	350623.668431141,	276963.341758769,	376659.066908011,	280208.201644618,	299011.810307279,	347881.649914261,	304804.276148251,	354140.442159242,	297577.470476564,	369013.453444279,	342425.115001835,	344230.596080637,	352928.701873343,	307616.284257507,	298634.310764378,	278128.01892035,	287575.078022483,	355882.380454943,	322944.076189626,	363781.805971695,	398048.237932872,	310256.468411648,	341124.312203455,	311328.460993228,	302880.616489425,	296513.384314058,	304445.37672747,	308823.524507046,	318609.641938095),
  DNA = c(2128800.34948927,	2161792.03497542,	2137106.03848348,	2080966.24118692,	1997683.41671565,	2002780.31225133,	1988922.90836014,	1914991.34151787,	2062256.61048089,	2068074.75438787,	2065842.75806008,	2066189.33561765,	2098999.30866258,	2053992.51463535,	2056188.52004091,	2081248.47961459,	2072758.14802907,	2049286.87801318,	2070404.44781516,	2082134.90006318,	2134881.06269515,	2141093.50600476,	2113932.09118523,	2126807.27629904,	2104528.66637072,	2070796.37762262,	2113910.47861006,	2084861.04940761,	2055734.90980933,	2088315.43659103,	2094359.44399986,	2073500.70388984,	2077115.60250976,	2061794.7010326,	2081982.84586638,	2057613.62977519,	2374160.35473021,	2142965.45896832,	2110611.80111109,	2097204.75258612,	1988806.19590372,	1940044.92597392,	2001853.82734762,	1959139.29867794,	2065615.74775278,	2059408.82645089,	2056624.60037196,	2063086.25144923,	2072933.46784475,	2104483.52935856,	2057095.56868244,	2067039.2743085,	2095673.79214113,	2077711.57351031,	2085082.74602288,	2100999.59845964,	2059583.84960386,	2071815.01705802,	2049988.49774249,	2058731.67605747,	2062719.09555333,	2065507.90733339,	2074692.05604424,	2047376.58085985,	2079353.51399334,	2102539.70122282,	2132336.2408236,	2101873.07112957,	2106918.53845382,	2106527.44028428,	2041452.01978421,	2077700.57096506,	2073918.42090824,	2065602.89089082,	2062348.87565966,	2069066.59177609,	2093284.95779904,	2062061.63161589,	2077924.30301342,	2091164.94785627),
  SATELLITE = c(112078.218171013,	112261.127881991,	109822.653395305,	108319.011290447,	120938.567331701,	124614.353474342,	123553.323943806,	121457.929807619,	121011.420265509,	122083.440200248,	118596.365414173,	126197.860989848,	114659.269110932,	115652.679257273,	117887.4889797,	116264.087295487,	118037.536136759,	116966.65686469,	116401.932894693,	120736.146934793,	121856.609173608,	115560.66768733,	116136.462931935,	116627.866510078,	113967.947930023,	115157.302494374,	115626.679437972,	112847.022338116,	123822.303554068,	117584.239385159,	113919.292366408,	120663.935405361,	115075.595108781,	116803.141784668,	114257.665673354,	119357.662099962,	158721.449077624,	122571.127415022,	111696.449638246,	117296.246940187,	114831.820664305,	124091.405662691,	119396.050074681,	116155.915367505,	127335.307025999,	124228.514297119,	124883.505437649,	122038.533798561,	115508.197369873,	121158.693720719,	117378.669222455,	116989.706969702,	118253.200519467,	122114.592330201,	124697.986038451,	119687.27985586,	116004.450402571,	122553.284795495,	115636.288470371,	113504.686350157,	117035.078585694,	121303.52633059,	117165.659762965,	122023.107960173,	113554.481867581,	108045.241904751,	116661.153273314,	114158.640138514,	107727.294522268,	115743.979540052,	115500.13156686,	116089.192819991,	114043.117289783,	117403.931665194,	116616.852366405,	115541.544461697,	118719.195699928,	120732.019580937,	115209.467638405,	119487.621780493),
  SINE = c(31377.272672055,	30132.7891361424,	29790.3726313403,	29358.1965060814,	40408.2545421847,	39859.65092241,	38730.7919323771,	43288.6072627528,	34244.3818588339,	33273.646344081,	34132.9016197188,	33126.1604970074,	31927.1949589471,	32381.4953306635,	33912.5441778959,	34665.1449202147,	29106.009118309,	29713.0902045069,	30864.9116155925,	31038.9314491018,	30131.7638190899,	30260.8538722906,	30319.7005520439,	28710.5827474677,	30660.7758661474,	30418.771083416,	29766.5863541231,	29822.062144439,	32070.5893435523,	30231.5118835594,	30593.7054209905,	30665.7874158144,	30082.29461937,	34596.14988353,	32135.5646840887,	31578.0549363139,	41446.9959995512,	34858.3529714664,	32636.4508023559,	31039.8119426606,	37765.5728460116,	43699.1217278204,	39441.6401884613,	36067.8887862087,	32457.3242954421,	35635.6254927885,	35023.5811071167,	32918.6344899252,	31659.4326706297,	33293.3620928082,	33076.3343934216,	30640.3631911463,	30254.7841828176,	31699.7877994612,	29533.8517417534,	29904.5506524809,	31707.4928839813,	34489.6159060525,	33524.0413170978,	29716.9061232749,	34810.4525738952,	33718.2298157538,	32012.6838465009,	33041.3497856505,	29202.6081279729,	30623.93913951,	29098.4633588707,	28588.1794597829,	30330.3442061969,	30560.7498996117,	31584.9874349504,	29596.9919350805,	31031.7329195725,	31330.3331265398,	33527.2495880913,	32371.3231756449,	30868.1951341776,	29763.7750714081,	33465.6727025237,	31765.0040389435),
  RC = c(65163.767485776,	66296.3571170775,	66069.5086904498,	66551.5047401099,	70870.1894790459,	71615.7546519351,	68431.2935516944,	70169.4951540686,	69102.3181629871,	66253.7016910082,	68044.9518849354,	68530.1567588615,	66687.6677601832,	68541.7425697181,	69434.0026257892,	68799.2856079294,	68956.4945398541,	68585.0975221352,	68542.6347729953,	69353.4911683249,	71004.5663600374,	70515.0177852271,	66978.2693484886,	67996.062484335,	70088.4177384189,	67998.7523564959,	68790.4751081596,	69111.911516055,	73689.0503812343,	70001.8350927709,	67151.4860580805,	68186.4226013712,	73214.0305535257,	72487.9433703464,	70646.1843527907,	73709.7974270855,	87823.7021867134,	67540.7567948746,	63681.7652787385,	69191.0048781171,	65856.9114868712,	67915.2134675882,	68253.9772822126,	67404.4464741765,	69809.3740171165,	71829.0200421666,	68503.2793408382,	68250.3476972768,	65094.6954727052,	68291.9562823282,	69908.0401017005,	69790.9375802206,	68235.1708730561,	71117.9248462557,	68871.5947095983,	71208.5316953775,	71757.6936284677,	72149.7687538329,	72736.4751972551,	70036.4812745246,	73921.8817966749,	73451.9918275525,	73046.915927408,	75471.4984174966,	67831.3991491695,	66152.6893170406,	70842.9158809714,	69693.4803447198,	66829.8944425812,	69718.3580003668,	75179.6073922849,	71076.6374560374,	66913.2527970981,	72218.5578057499,	68668.5321488357,	71077.1703476654,	70861.4105782479,	72665.2688603354,	71442.5121960958,	72161.2011151279)
)

# Reshape the data from wide to long format
data_long <- gather(data, key = Variable, value = Value, -Treatment, -gene)

# Create a box plot with data points and facets by variable
p <- ggplot(data_long, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_boxplot() +
  geom_point(size = 1, alpha = 0.5) +
  facet_grid(. ~ Variable, scales = "free_y") +
  ggtitle("RNASeq TEtranscripts: Social stress TE class counts") +
  xlab("Social stress treatment") +
  ylab("Normalized Counts (Log10 Scale)") +
  scale_fill_manual(values = c("Pink", "Purple")) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold", color = "black"),
        axis.title = element_text(size = 16, face = "bold", color = "black"),
        axis.line = element_line(),
        axis.ticks = element_line()) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20)) +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size = 14, face = "bold"))


p


p + scale_y_log10(limits = c(3e+04, 1e+65)) # just to give a bit more space to manually add significance stars on the plot later


# mann whitney u
#spss kolmogorov smirnov used for tests of normality- data not normal distributed


data <- data.frame(
  Treatment = c("High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"High",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low", "Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low",	"Low"),
  
  #LINE = c(172051.010343351,	181798.171233951,	167239.777767372,	192813.461754451,	234507.912857316,	198465.146952479,	219371.804003757,	254715.971981223,	177872.599878578,	195369.486516681,	186575.549064863,	201324.637555469,	179216.10968571,	186421.800763887,	185120.128684641,	198193.776832608,	193190.522875243,	193527.751991868,	181011.622541169,	184507.011935077,	226235.773212631,	197894.210778839,	189854.351853873,	184072.883442402,	174914.796878524,	172060.599557008,	171793.431070693,	168941.949434443,	180536.392262371,	169744.986519629,	174388.797247282,	186705.018684362,	177380.640450672,	185643.461091638,	173591.46622382,	182114.696361727,	231930.833119505,	219266.959940662,	205004.876206464,	183653.148517671,	202136.890583854,	207700.262207727,	195273.538685284,	189363.016369216,	196400.954065268,	186033.660168918,	184038.385310046,	186215.735106491,	184769.094221855,	183421.014720214,	184304.639842552,	180285.926030174,	194996.863373054,	185060.206594937,	232103.284254028,	213714.890062587,	201696.466066612,	215539.518364105,	192318.224264503,	191506.288687014,	208302.467984097,	196633.330097499,	189652.449112463,	198504.001771491,	214790.24181259,	187171.442263759,	179995.988899247,	189170.561691902,	175121.940174235,	167420.623928208,	180161.10058392,	197874.731506271,	184892.33691378,	177169.710953803,	181255.393670272,	184596.950712525,	182702.720566856,	178633.241812634,	182486.315808998,	174698.646682999)
  #LTR = c(277780.208508856,	271662.579920513,	296765.936026165,	284761.0094366,	377902.905743328,	381430.290894118,	380749.370267615,	489972.54731661,	331035.299948178,	308977.340152634,	377095.568277416,	324245.090479094,	292806.353632439,	343807.728534855,	356311.348923484,	392302.309423938,	321100.436413208,	319298.147192926,	320190.350886054,	290906.780137328,	297408.829586074,	255594.055233358,	309829.121959588,	296220.638665861,	294594.997037001,	310879.856844734,	316555.937358353,	318786.886085994,	326885.371625091,	319104.167137172,	301124.746096668,	295943.277073925,	327896.431805257,	338734.92982271,	297646.111436376,	303651.843397811,	436543.626847255,	321077.594390848,	309852.342250699,	317763.019481653,	360952.993110389,	413393.838420643,	362947.175144834,	386651.88078083,	345211.088563561,	335475.845214327,	319404.718973342,	332128.615494627,	296204.415940727,	319459.834446284,	335160.67234884,	350623.668431141,	276963.341758769,	376659.066908011,	280208.201644618,	299011.810307279,	347881.649914261,	304804.276148251,	354140.442159242,	297577.470476564,	369013.453444279,	342425.115001835,	344230.596080637,	352928.701873343,	307616.284257507,	298634.310764378,	278128.01892035,	287575.078022483,	355882.380454943,	322944.076189626,	363781.805971695,	398048.237932872,	310256.468411648,	341124.312203455,	311328.460993228,	302880.616489425,	296513.384314058,	304445.37672747,	308823.524507046,	318609.641938095)
  #DNA = c(2128800.34948927,	2161792.03497542,	2137106.03848348,	2080966.24118692,	1997683.41671565,	2002780.31225133,	1988922.90836014,	1914991.34151787,	2062256.61048089,	2068074.75438787,	2065842.75806008,	2066189.33561765,	2098999.30866258,	2053992.51463535,	2056188.52004091,	2081248.47961459,	2072758.14802907,	2049286.87801318,	2070404.44781516,	2082134.90006318,	2134881.06269515,	2141093.50600476,	2113932.09118523,	2126807.27629904,	2104528.66637072,	2070796.37762262,	2113910.47861006,	2084861.04940761,	2055734.90980933,	2088315.43659103,	2094359.44399986,	2073500.70388984,	2077115.60250976,	2061794.7010326,	2081982.84586638,	2057613.62977519,	2374160.35473021,	2142965.45896832,	2110611.80111109,	2097204.75258612,	1988806.19590372,	1940044.92597392,	2001853.82734762,	1959139.29867794,	2065615.74775278,	2059408.82645089,	2056624.60037196,	2063086.25144923,	2072933.46784475,	2104483.52935856,	2057095.56868244,	2067039.2743085,	2095673.79214113,	2077711.57351031,	2085082.74602288,	2100999.59845964,	2059583.84960386,	2071815.01705802,	2049988.49774249,	2058731.67605747,	2062719.09555333,	2065507.90733339,	2074692.05604424,	2047376.58085985,	2079353.51399334,	2102539.70122282,	2132336.2408236,	2101873.07112957,	2106918.53845382,	2106527.44028428,	2041452.01978421,	2077700.57096506,	2073918.42090824,	2065602.89089082,	2062348.87565966,	2069066.59177609,	2093284.95779904,	2062061.63161589,	2077924.30301342,	2091164.94785627)
  #SATELLITE = c(112078.218171013,	112261.127881991,	109822.653395305,	108319.011290447,	120938.567331701,	124614.353474342,	123553.323943806,	121457.929807619,	121011.420265509,	122083.440200248,	118596.365414173,	126197.860989848,	114659.269110932,	115652.679257273,	117887.4889797,	116264.087295487,	118037.536136759,	116966.65686469,	116401.932894693,	120736.146934793,	121856.609173608,	115560.66768733,	116136.462931935,	116627.866510078,	113967.947930023,	115157.302494374,	115626.679437972,	112847.022338116,	123822.303554068,	117584.239385159,	113919.292366408,	120663.935405361,	115075.595108781,	116803.141784668,	114257.665673354,	119357.662099962,	158721.449077624,	122571.127415022,	111696.449638246,	117296.246940187,	114831.820664305,	124091.405662691,	119396.050074681,	116155.915367505,	127335.307025999,	124228.514297119,	124883.505437649,	122038.533798561,	115508.197369873,	121158.693720719,	117378.669222455,	116989.706969702,	118253.200519467,	122114.592330201,	124697.986038451,	119687.27985586,	116004.450402571,	122553.284795495,	115636.288470371,	113504.686350157,	117035.078585694,	121303.52633059,	117165.659762965,	122023.107960173,	113554.481867581,	108045.241904751,	116661.153273314,	114158.640138514,	107727.294522268,	115743.979540052,	115500.13156686,	116089.192819991,	114043.117289783,	117403.931665194,	116616.852366405,	115541.544461697,	118719.195699928,	120732.019580937,	115209.467638405,	119487.621780493)
  #SINE = c(31377.272672055,	30132.7891361424,	29790.3726313403,	29358.1965060814,	40408.2545421847,	39859.65092241,	38730.7919323771,	43288.6072627528,	34244.3818588339,	33273.646344081,	34132.9016197188,	33126.1604970074,	31927.1949589471,	32381.4953306635,	33912.5441778959,	34665.1449202147,	29106.009118309,	29713.0902045069,	30864.9116155925,	31038.9314491018,	30131.7638190899,	30260.8538722906,	30319.7005520439,	28710.5827474677,	30660.7758661474,	30418.771083416,	29766.5863541231,	29822.062144439,	32070.5893435523,	30231.5118835594,	30593.7054209905,	30665.7874158144,	30082.29461937,	34596.14988353,	32135.5646840887,	31578.0549363139,	41446.9959995512,	34858.3529714664,	32636.4508023559,	31039.8119426606,	37765.5728460116,	43699.1217278204,	39441.6401884613,	36067.8887862087,	32457.3242954421,	35635.6254927885,	35023.5811071167,	32918.6344899252,	31659.4326706297,	33293.3620928082,	33076.3343934216,	30640.3631911463,	30254.7841828176,	31699.7877994612,	29533.8517417534,	29904.5506524809,	31707.4928839813,	34489.6159060525,	33524.0413170978,	29716.9061232749,	34810.4525738952,	33718.2298157538,	32012.6838465009,	33041.3497856505,	29202.6081279729,	30623.93913951,	29098.4633588707,	28588.1794597829,	30330.3442061969,	30560.7498996117,	31584.9874349504,	29596.9919350805,	31031.7329195725,	31330.3331265398,	33527.2495880913,	32371.3231756449,	30868.1951341776,	29763.7750714081,	33465.6727025237,	31765.0040389435)
  RC = c(65163.767485776,	66296.3571170775,	66069.5086904498,	66551.5047401099,	70870.1894790459,	71615.7546519351,	68431.2935516944,	70169.4951540686,	69102.3181629871,	66253.7016910082,	68044.9518849354,	68530.1567588615,	66687.6677601832,	68541.7425697181,	69434.0026257892,	68799.2856079294,	68956.4945398541,	68585.0975221352,	68542.6347729953,	69353.4911683249,	71004.5663600374,	70515.0177852271,	66978.2693484886,	67996.062484335,	70088.4177384189,	67998.7523564959,	68790.4751081596,	69111.911516055,	73689.0503812343,	70001.8350927709,	67151.4860580805,	68186.4226013712,	73214.0305535257,	72487.9433703464,	70646.1843527907,	73709.7974270855,	87823.7021867134,	67540.7567948746,	63681.7652787385,	69191.0048781171,	65856.9114868712,	67915.2134675882,	68253.9772822126,	67404.4464741765,	69809.3740171165,	71829.0200421666,	68503.2793408382,	68250.3476972768,	65094.6954727052,	68291.9562823282,	69908.0401017005,	69790.9375802206,	68235.1708730561,	71117.9248462557,	68871.5947095983,	71208.5316953775,	71757.6936284677,	72149.7687538329,	72736.4751972551,	70036.4812745246,	73921.8817966749,	73451.9918275525,	73046.915927408,	75471.4984174966,	67831.3991491695,	66152.6893170406,	70842.9158809714,	69693.4803447198,	66829.8944425812,	69718.3580003668,	75179.6073922849,	71076.6374560374,	66913.2527970981,	72218.5578057499,	68668.5321488357,	71077.1703476654,	70861.4105782479,	72665.2688603354,	71442.5121960958,	72161.2011151279)
)

# Separate the LINE values for each treatment group
rc_high <- data$RC[data$Treatment == "High"]
rc_low <- data$RC[data$Treatment == "Low"]

# Perform the Mann-Whitney U test
result <- wilcox.test(rc_high, rc_low, alternative = "two.sided")
print(result)


# box plot with mann whitney u testing
# Load the necessary libraries
library(ggplot2)

# Create a data frame from the image data
data <- data.frame(
  Treatment = c("High", "High", "High", "High", "High", "High", "High", "High", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "High", "High", "High", "High", "High", "High", "High", "High", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "High", "High", "High", "High", "Low", "Low", "Low", "Low", "High", "High", "High", "High", "High", "High", "High", "High", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "High", "High", "High", "High", "High", "High", "High", "High"),
  gene = c("E100-1", "E100-2", "E100-3", "E100-4", "E101-1", "E101-2", "E101-3", "E101-4", "E102-1", "E102-2", "E102-3", "E102-4", "E103-1", "E103-2", "E103-3", "E103-4", "E104-1", "E104-2", "E104-3", "E104-4", "E105-1", "E105-2", "E105-3", "E105-4", "E106-1", "E106-2", "E106-3", "E106-4", "E107-1", "E107-2", "E107-3", "E107-4", "E108-1", "E108-2", "E108-3", "E108-4", "E109-1", "E109-2", "E109-3", "E109-4", "E110-1", "E110-2", "E110-3", "E110-4", "E111-1", "E111-2", "E111-3", "E111-4", "E112-4", "E112-5", "E112-7", "E112-8", "E113-10", "E113-11", "E113-13", "E113-9", "E114-1", "E114-3", "E114-4", "E114-5", "E115-1", "E115-2", "E115-4", "E115-5", "E80-1", "E80-2", "E80-3", "E80-4", "E81-1", "E81-2", "E81-4", "E81-5", "E82-1", "E82-2", "E82-3", "E82-5", "E83-1", "E83-3", "E83-4", "E83-5"),
  value = c(421.172611350585,	487.167670393007,	333.425331389698,	224.200416315901,	689.71795130167,	564.155539947648,	466.104780320239,	543.423683028866,	1649.07516628642,	800.27148806119,	1268.37624379016,	242.409126132341,	537.834835264754,	717.822008375384,	494.146724486718,	515.701061974672,	1631.43695964258,	1523.73273032068,	1222.32214886215,	2611.38263052907,	492.616092868772,	533.303158063629,	472.271125030913,	541.935039890007,	1419.50857645487,	1639.87119045591,	1247.23411774184,	1345.06502379872,	625.165963080654,	1224.17774646106,	424.882523712694,	440.494625215946,	397.685784665952,	478.216629706535,	525.332435488413,	594.433597640564,	218.907244529264,	269.059034997562,	195.64275079184,	381.917500883641,	121.743891427299,	113.95389844184,	188.105157612644,	352.955683945291,	276.577898718085,	318.909848155069,	160.570757880699,	320.590921383672,	449.965139434238,	1248.62934897101,	985.539075547919,	726.97783831206,	1343.81383529391,	841.409039079279,	813.934636208223,	1016.30863232746,	346.06442490521,	142.436761804881,	141.023998250619,	134.755305968953,	502.863040739236,	409.286654604506,	350.540034293381,	273.942355762314,	2510.81931261142,	2260.43746901733,	2539.29150951435,	2903.51448961465,	1896.27086718087,	2164.71986478154,	1838.51780017382,	1808.97945713209,	2365.82981589229,	2911.43315507587,	2918.21831664749,	2369.77469786511,	1590.0307954493,	1400.12185294625,	1363.07705268363,	1585.33140700549)
)



# Remove non-finite values from the data
data <- data[data$value > 0, ]


# Create a box plot with data points
p <- ggplot(data, aes(x = Treatment, y = value, fill = Treatment)) +
  geom_boxplot() +
  geom_point(size = 1, alpha=0.5) +
  ggtitle("Box Plot of Gene Expression by Treatment (Log10 Scale)") +
  xlab("Treatment") +
  ylab("Normalized Counts (Log10 Scale)") +
  scale_fill_manual(values = c("Pink", "Purple")) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, face="bold", color="black"),
        axis.text.y = element_text(size = 14, face="bold", color="black"),
        axis.title = element_text(size = 14, face="bold", color="black"),
        #legend.title = element_text(size = 14, face="bold", color="black"),
        #legend.text = element_text(size=14, face="bold", color="black"),
        #panel.grid.major = element_blank(), # Remove grid lines
        #panel.grid.minor = element_blank(),
        axis.line = element_line(), # Add axis lines
        axis.ticks = element_line()) +   # Add tick marks 
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  theme(legend.position = "none")

p
p + ggtitle("RNA-seq normalized gene counts for: Metazoa SRP") +
  scale_y_log10(limits = c(100, 4500))
p


p + scale_y_log10(limits = c(1e+03, 1e+06)) +# just to give a bit more space to manually add significance stars on the plot later
  
  # mann whitney u test
  # Separate the LINE values for each treatment group
  gene_high <- data$value[data$Treatment == "High"]
gene_low <- data$value[data$Treatment == "Low"]

# Perform the Mann-Whitney U test
result <- wilcox.test(gene_high, gene_low, alternative = "two.sided")
print(result)

#heatmap with pheatmap package

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

data <- read.csv("sig_de_mirnas_40_volcano.csv")#, comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names

head(mat_data)
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
         fontsize = 16,
         fontface = "bold",
         color = colorRampPalette(inferno(10))(10000))



####Getting TE counts from TEtranscripts

data <- read.table("socstress_highvlow.cntTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",36),rep("CGroup",44)))
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)

keep <- rowSums(counts(dds)) >= 10

dds <- dds [keep,]

dds
dds <- DESeq (dds)
resultsNames(dds)
res <- results(dds)


write.table(res, file="socstress_highvlow_gene_TE_analysis_10.txt", sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 0.000000)), ]
write.table(resSig, file="socstress_highvlow_sigdiff_gene_TE_10.txt",sep="\t", quote=F)
