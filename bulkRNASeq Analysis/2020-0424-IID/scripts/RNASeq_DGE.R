# ---
# title: "Basic RNASeq Analysis - DGE"
# author: "Yaoyu E. Wang"
# output: html_document
# ---
# 
# 
# 
# # Begin RNA-Seq DGE Analysis
# We first load the RNA-Seq count matrix from 'RNASeqData.Rdata'.  This is a simple 3 vs 3 RNA-Seq experiment with Control vs Experiment in count matrix by gene.    
# 
# ```{r load, echo = TRUE, include=TRUE}
# One can set up the path
# Set working directory to PROJECT_DIR
setwd("~/Dropbox (Harvard University)/qBRC-Share/Presentations/QBRC-Workshops/Intro to RNASeq Analysis/2020-0424-IID/RNASeq_Exploratory")
# setwd("~/modules/rnaseq_dge/")
load('RNASeq_exploratory.RData') 


ls()                       # list all current variables

head(count_data)           # The head(count_data) provides the first 6 row of the 'count_data'
head(groups)           # The head(count_data) provides the first 6 row of the 'count_data'


## Principle component analysis
# 
# Once the matrix is loaded,  we can perform principle component analysis to visualize sample distribution.  
# We first group samples into two vectors/variables: CONTROL and TREATED.
# 
# Calculate PCA use 'prcomp' command.  Since prcomp compute PCA by the rows, we will need to tranpose the count matrix 
# such that the samples are represented by rows and genes by the columns.
# 
# ```{r pca, echo=TRUE}
CONTROL=groups$control$sample
TREATED=groups$treated$sample

pca=prcomp(t(count_data))        # compute PCA on transposed count_data
plot(pca$x[,1], pca$x[,2], xlab='pc1', ylab='pc2')
points(pca$x[CONTROL,1], pca$x[CONTROL,2],
       col="red", pch=18)
points(pca$x[TREATED,1], pca$x[TREATED,2],
       col="blue", pch=18)

# # ## Hierarchical Clustering
# # 
# # We first calculate similarity matrix using Euclidean distance.  The matrix is computed by 'dist' function. 'hclust' 
# # function computes hierarhical cluster tree and save it into an variable object 'hctree' that can be ploted by using generic 'plot' function.
# 
# ```{r hclust, echo=TRUE}

# compute similarity matrix on the transposed count_data matrix using Euclidean distance
dist_matrix=dist(t(count_data), method="euclidean")
# compute hierarchical clustering tree and save it into variable called 'hctree'
hctree=hclust(dist_matrix)
# plot hctree
plot(hctree)

# 
# ## Differential Gene Expression
# 
# We will use DESeq2 to perform differential gene expression.   Since DESeq2 is a specialized bioconductor package/library, 
# we will need to install it before loading in the package to use.
# 
# The following command check if DESeq2 and Bioconductor Manager packages have been installed, and install the packages if they have not been installed.
# ```{r DESeq2, include=TRUE}
#
# Check for Bioconductor, install it if it is not yet installed
#
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

# Check for DESeq2, install it if it is not yet installed
if(!requireNamespace('DESeq2', quietly = TRUE))
    BiocManager::install('DESeq2', update=FALSE)

library(DESeq2)

# We can then use the functions within DESeq2 to perform DGE analysis
# First, we define which sampels are control group and which samples are treated group
# The samples are ordered as YEW1_Control, YEW2_Control, YEW3_Control, YEW1_Treated, 
# YEW2_Treated, YEW3_Treated, so we can use the following to define conditions of 
# each sample
# 
# ```{r condition, echo=TRUE}
#

condition <- factor(c("control", "control", "control",
                      "treated", "treated", "treated"))
# alternatively, since the count_data is ordered in the same way as groups variable
# we can declare condition from groups:
# 
# condition <- factor(c(as.character(groups$control$group), as.character(groups$treated$group)))
#


# We then use the 'condition' to perform contrast on the count_data.  Noted that count_data is still raw count matrix,  D
# ESeq2 performs normalization and DGE together in one function.
# 
# If you want to know how to use a specific function in R,  just type '?function' in Console.
# For example type:  
# 
# > ?DESeqDataSetFromMatrix
#
# ```{r DGE, include=FALSE}
#
# This function generates model and set up contrast using linear model

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = DataFrame(condition),
                              design = ~condition)

# This function runs the model and save all model parameters and estimation into dds object
dds <- DESeq(dds)

# Returns DGE gene list results in data frame format
dge_results=results(dds, contrast=c("condition", "treated", "control"))

# make sure there is no missing values
dge_results=na.omit(dge_results)  

# filter for results with padj<0.001 and absolute log2FoldChange>2
filtered_results=subset(dge_results, padj<0.001 & abs(log2FoldChange)>2)

# Obtain normalized count values
norm_count <- counts(dds, normalized = T)

# output filtered_results
filtered_results

#
# We can compute PCA on normalized matrix again
#```{r norm_pca, echo=TRUE}
#
CONTROL=groups$control$sample
TREATED=groups$treated$sample

pca=prcomp(t(norm_count))        # compute PCA on transposed normalized_data
plot(pca$x[,1], pca$x[,2], xlab='pc1', ylab='pc2')
points(pca$x[CONTROL,1], pca$x[CONTROL,2],
       col="red", pch=18)
points(pca$x[TREATED,1], pca$x[TREATED,2],
       col="blue", pch=18)

# compute similarity matrix on the transposed count_data matrix using Euclidean distance
dist_matrix=dist(t(norm_count), method="euclidean")
# compute hierarchical clustering tree and save it into variable called 'hctree'
hctree=hclust(dist_matrix)
# plot hctree
plot(hctree)


# 
# ## Heatmap
# 
# Generate heatmap using filtered dge results 'filtered_results'.  The heatmap function is core function, but we want to have better coloring, so we install and load gplots. 
# 
# ```{r heatmap, include=TRUE}
# Check and install "gplots" if not installed

if(!requireNamespace('gplots', quietly = TRUE))
  install.packages('gplots')
library(gplots)

DGE_matrix=norm_count[rownames(filtered_results),]
heatmap(DGE_matrix,col=colorpanel(50, 'blue', 'white', 'red'),
        Rowv=TRUE, Colv=TRUE)

# call the commmand again to write it into a file.
png("heatmap.png")
heatmap(DGE_matrix,col=colorpanel(50, 'blue', 'white', 'red'),
        Rowv=TRUE, Colv=TRUE)
dev.off()


## Generate Volcano Plot
#
# Again, we first check for the package 'EnhancedVolcano', and install it if it is not already installed. 
#
# ```{r volcano, echo=FALSE}
if (!requireNamespace('EnhancedVolcano', quietly = TRUE))
    BiocManager::install('EnhancedVolcano', update=FALSE)


#Load the library and run EnhancedVolcano function to generate volcano plot.
library(EnhancedVolcano)

volcano_plot<-EnhancedVolcano(dge_results,
    lab = rownames(dge_results),
    x = 'log2FoldChange',
    y = 'padj',
    xlim = c(-5, 8),
    ylim = c(0, 20),
    #pCutoff = 10e-16,
    FCcutoff = 1.5
)

volcano_plot

 
# ## Gene Set Erichment Analysis
# 
# There are many different packages retreiving MSigDB data to perform GSEA analysis. Most of the packages perform these tasks in very similar ways with different application programming interface (API).  We will use the following packages:
# 
# - migdbr (https://github.com/igordot/msigdbr)
# - fgsea  (https://bioconductor.org/packages/release/bioc/html/fgsea.html)
# 
# 
# First, we install **misgdbr**
# ```{r misgdbr_install, include=TRUE}
# install misgdbr from CRAN
if(!require(msigdbr)) install.packages("msigdbr")
library(msigdbr)
# 
# The **misgdbr** retrieve data from MSigDB database hosted at Broad Institute (https://www.gsea-msigdb.org/gsea/msigdb/index.jsp).
# The function **misgdbr** runs as the following:
# 
# ####Usage
# msigdbr(species = "Homo sapiens", category = NULL, subcategory = NULL)
# 
# ####Arguments
# |Parameters|Description|
# |:---:|---|
# |**species**|species name, such as Homo sapiens, Mus musculus, etc.|
# |**category**|collection, such as H, C1, C2, C3, C4, C5, C6, C7.|
# |**subcategory**|sub-collection, such as CGP, MIR, BP, etc.|
# 
# ```{r gsea, include=TRUE}

m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
# reformat the data frame into list for fgsea input
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name) 


# install fgsea
if(!require(fgsea)) BiocManager::install("fgsea")
library(fgsea)

# reformat
gene_stat=dge_results$stat
names(gene_stat)=rownames(dge_results)
# sort genes
gene_stat=gene_stat[order(gene_stat)]

fgsea_results=fgsea(pathways=m_list, gene_stat, nperm=50, minSize=20)
fgsea_results=fgsea_results[order(fgsea_results$pval),]
head(fgsea_results)

#
# Now write out the GSEA results
# 
# ```{r output}
library(writexl)    # load  library to write excel file

# Current directory is Data directory, change to output
if(!dir.exists('outputs')){
  dir.create('outputs')
}

# using ggplots's ggsave function to save the plot
print("Save the volcano plot in ggplot format")
ggsave("outputs/volcano.png", plot=volcano_plot)


print('Write the final DGE data to Excel file.')
# start a new excel file and write final results 
write_xlsx(as.data.frame(dge_results), "outputs/dge_results.xlsx", col_names=TRUE)


print('Write the final GSEA data to Excel file.')
# start a new excel file and write final results 
write_xlsx(fgsea_results, "outputs/gsea_results.xlsx", col_names=TRUE)
