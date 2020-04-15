# Load library
if(!requireNamespace('gplots', quietly = TRUE))
  install.packages('gplots')
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
if(!requireNamespace('DESeq2', quietly = TRUE))
  BiocManager::install('DESeq2', update=FALSE)
if (!requireNamespace('EnhancedVolcano', quietly = TRUE))
  BiocManager::install('EnhancedVolcano', update=FALSE)

library(EnhancedVolcano)
library(gplots)
library(DESeq2)



# Begin RNA-Seq Analysis

#We first load the RNA-Seq count matrix from 'RNASeqData.Rdata'.  This is a simple 3 vs 3 RNA-Seq experiment with Control vs Experiment in count matrix by gene.    

load('RNASeqData.Rdata')   # Load data set which contains a variable called 'count_data'
head(count_data)           # The head(count_data) provides the first 6 row of the 'count_data'


## Principle component analysis
# Once the matrix is loaded,  we can perform principle component analysis to visualize sample distribution.  
# We first group samples into two vectors/variables: CONTROL and TREATED

CONTROL=c("YEW1_Control", "YEW2_Control", "YEW3_Control")
TREATED=c("YEW4_Treated", "YEW5_Treated", "YEW6_Treated")

#Calculate PCA use 'prcomp' command.  Since prcomp compute PCA by the rows, we will need to tranpose the count matrix such that the samples are represented by rows and genes by the columns.

pca=prcomp(t(count_data))        # compute PCA on transposed count_data

#Plot first two dimension of pca results, e.g. pc1 vs pc2 and color CONTROL samples in red and treated samples in blue.
plot(pca$x[,1], pca$x[,2], xlab='pc1', ylab='pc2')
points(pca$x[CONTROL,1], pca$x[CONTROL,2],
       col="red", pch=18)
points(pca$x[TREATED,1], pca$x[TREATED,2],
       col="blue", pch=18)


## Hierarchical Clustering

# We first calculate similarity matrix using Euclidean distance.  The matrix is computed by 'dist' function. 'hclust' function computes hierarhical cluster tree and save it into an variable object 'hctree' that can be ploted by using generic 'plot' function.

# compute similarity matrix on the transposed count_data matrix using Euclidean distance
dist_matrix=dist(t(count_data), method="euclidean")
# compute hierarchical clustering tree and save it into variable called 'hctree'
hctree=hclust(dist_matrix)
# plot hctree
plot(hctree)

## Differential Gene Expression

# We will use DESeq2 to perform differential gene expression.   Since DESeq2 is a specialized bioconductor package/library, we will need to install it before loading in the package to use.

# The following command check if DESeq2 and Bioconductor Manager packages have been installed, and install the packages if they have not been installed.


#Once the packages are installed, we can call 'library' function to load 'DESeq2' library to use its functions.


# We can then use the functions within DESeq2 to perform DGE analysis
# First, we define which sampels are control group and which samples are treated group
# The samples are ordered as YEW1_Control, YEW2_Control, YEW3_Control, YEW1_Treated, 
# YEW2_Treated, YEW3_Treated, so we can use the following to define conditions of 
# each sample

condition <- factor(c("control", "control", "control",
                      "treated", "treated", "treated"))

# We then use the 'condition' to perform contrast on the count_data.  Noted that count_data is still raw count matrix,  DESeq2 performs normalization and DGE together in one function.
# If you want to know how to use a specific function in R,  just type '?function' in Console.
# For example type:  
  
#  > ?DESeqDataSetFromMatrix


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

## Heatmap
# Generate heatmap using filtered dge results 'filtered_results'.  The heatmap function is core function such that no special package is needed.

DGE_matrix=norm_count[rownames(filtered_results),]
heatmap(DGE_matrix,col=colorpanel(50, 'blue', 'white', 'red'),
        Rowv=TRUE, Colv=TRUE)

## Generate Volcano Plot

# Again, we first check for the package 'EnhancedVolcano', and install it if it is not already installed. 

a=EnhancedVolcano(dge_results,
                  lab = rownames(dge_results),
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlim = c(-5, 8),
                  ylim = c(0, 20),
                  FCcutoff = 1.5
)

a

