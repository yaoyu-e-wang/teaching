#
# This script is the final result of all the analyses we will perform on the Data Wrangling with R
# real data section
#

# First, set up the proper directory path and load all the data:
setwd("~/modules/data_wrangling")

print('Loading data...')
dge_results <- read.table('data/differential_results.csv', sep=',', header=T)
expressions = read.table('data/gene_expression.tsv', header=T)
annotations = read.delim('data/gene_annotations.tsv')
pathways = read.table('data/my_pathway_genes.txt', sep='\t', col.names=c('gene_name','pathway'))
mutations = read.delim('data/mutations.tsv')
print('Done loading data.')

###########################
# Now we have all the data store as variables
#
print("Show the first six rows of each variable")
head(dge_results)
head(expressions)
head(annotations)
head(pathways)
head(mutations)
print("Done showing the first row of each data")

########################################
# The data is like excel worksheet. Let's look at the data
########################################
head(expressions)
# These two lines return the same number
expressions[2,2]
expressions[2, 'SW1_Control']

# Select First Column:
expressions[1]  
expressions['gene']  
expressions[,1]  
expressions$gene

# select First Row:
expressions[1,]

# select first three columns
expressions[1:3]

# select all columns except first:
expressions[2:ncol(expressions)]



#########################################################################################
# Boolean Comparison: AND, OR, NOT, Equal
#
# select gene annotations that are on chromosome 7
df=read.table("data/demo_annotations.tsv", header=T)
chroms=df$chrom
is_chr7 = chroms == 'chr7'
is_chr7
df[is_chr7,]

is_chr3 <- chroms == 'chr3'
df[is_chr3,]

# select gene that are annotated as oncogenes use %in%
oncogenes = c('KRAS', 'TP53')  
is_oncogene = df$name %in% oncogenes
df[is_oncogene,]

# Combine multiple criteria with AND (&) operation
selection_criteria = is_oncogene & is_chr7
data_subset = df[selection_criteria, ]
data_subset


#########################################################################################
print('Filtering for ras signaling genes that are significantly upregulated.')
# Get a vector of gene names that are in the pathway of interest
ras_genes = pathways[pathways['pathway'] == 'ras_signaling', 'gene_name']

# Create a boolean (True or False) vector for those RAS genes
is_ras_gene = dge_results$gene %in% ras_genes

# Create a boolean vector for whether the gene is significantly changed
# at a p < 0.05 threshold
is_significant = dge_results$padj < 0.05

# Create a boolean vector for whether the gene is upregulated
is_upregulated = dge_results$log2FoldChange > 0

# Keep only rows that pass all 3 "tests" (are True for all 3)
selected_rows = dge_results[is_ras_gene & is_significant & is_upregulated, ]

# Just in case, remove missing data
selected_rows = na.omit(selected_rows)

# We now have the rows/genes-- now, keep only a subset of the columns:
upregulated_ras_genes = selected_rows[c('gene','baseMean','log2FoldChange','padj')]

print('Done filtering genes.')

#########################################################################################

print('Merge with gene coordinates.')
# Add in the coordinate information by merging
ras_up_genes_w_coords = merge(upregulated_ras_genes, annotations, by.x = 'gene', by.y='name')

print('Merge with gene expression')
ras_up_genes_w_coords_expression = merge(ras_up_genes_w_coords, expressions, by='gene')




#########################################################################################

# Select the first row from mutation and change chrom from '5' to 'chr5'
i=1
chrom =mutations[i,'chrom']
chrom_w_prefix = paste('chr', chrom, sep='')  
mutations[i,'chrom'] =chrom_w_prefix
                 

#########################################################################################

# Use flow control - the “for loop”
print('Example using the for loop flow control')
my_vector <- c(10,11,12,13) 
for (item in my_vector){
  print(item)
  #do actual operations here
}

########################################################################################

print('Changing chromosome names in mutation dataframe')
# We need to add the 'chr' prefix to the chromosome names that are in the mutations dataframe.
# This way, all the chromosome names are consistent.  
# The method below is a slow way to do this, but easier to understand.
for ( i in 1:nrow(mutations) )
{
	chrom = mutations[i,'chrom']
	chrom_w_prefix = paste('chr', chrom, sep='')
	mutations[i, 'chrom'] = chrom_w_prefix
}

#########################################################################################

print('Filter to keep only genes that are mutated.')

# Again, this a slow, but clear way to do this.
ras_up_mutated_genes = data.frame()					        # define an empty result variable
for ( i in 1:nrow(ras_up_genes_w_coords) )
{
  gene_info = ras_up_genes_w_coords[i,]       		  # get row i and save as gene_info
  same_chrom = gene_info$chrom == mutations$chrom		# get mutation row with same chrom
  past_start = mutations$pos >gene_info$start			  # get position > then gene_info start
  before_end = mutations$pos <gene_info$end			    # get position < then gene_info start	
  overlap = same_chrom & past_start & before_end		# get mutations fit all 3 conditions
  if ( any(overlap) )							                  # if there is any mutation quality
  {
    ras_up_mutated_genes = rbind(ras_up_mutated_genes, gene_info) # use rbind() to add to result variable
  }
}

#########################################################################################

print('Merge to add in the expression data')
# note no by.x or by.y since the 'gene' column is in each dataframe
final_data = merge(ras_up_mutated_genes, expressions)
final_data

#########################################################################################
# get current working directory. 
getwd()
# Current directory is Data directory, change to output
if(!dir.exists('outputs')){
  dir.create('outputs')
}
# write final data out as a tab-seperated-value file (.tsv)
print('Write the final data to file.')
write.table(final_data, 'outputs/wrangling_data.tsv', sep='\t', quote=F)

# if package not installed
# run: 
# install.packages("xlsx")
if(!require(readxl)) install.packages("readxl")
if(!require(writexl)) install.packages("writexl")

library(readxl)    # load  library to write excel file
library(writexl)    # load  library to write excel file



print('Write the final data to Excel file.')
# start a new excel file and write final results 
write_xlsx(final_data, "outputs/wrangling_results.xlsx", col_names=TRUE)
# add a new sheet with new data onto the file we just created
write_xlsx(ras_up_genes_w_coords, "outputs/ras_up_genes_w_coords.xlsx", col_names=TRUE)

#########################################################################################




