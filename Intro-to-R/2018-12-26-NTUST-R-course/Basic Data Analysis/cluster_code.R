library(reshape2)
library(tidyverse)

setwd('Basic Data Analysis')
load('data/IHME_GBD_2017_PreProcessed_DataSets.Rdata')

# head and tail
head(raw_dt)
head(matrix_dt)
head(perc_dt)
head(location_code)

# perform PCAs
pca=princomp(matrix_dt)  # compute PCA
png('results/pca.png', res=150)
plot(pca$scores[,1], pca$scores[,2],
     xlab='pc1', ylab='pc2')
points(pca$scores[regions$`South Asia`$country,1], 
       pca$score[regions$`South Asia`$country,2], 
       col="red", pch=18)
points(pca$scores[regions$`Western Europe`$country,1], 
       pca$score[regions$`Western Europe`$country,2], 
       col="green", pch=18)
points(pca$scores[regions$`Central Sub-Saharan Africa`$country,1], 
       pca$score[regions$`Central Sub-Saharan Africa`$country,2], 
       col="blue", pch=18)
dev.off()

png('results/hctree.png', width=1200, height=2400)
hctree=hclust(dist(matrix_dt, method='euclidean'))
par(mar=c(3,1,1,15)) 
plot(as.dendrogram(hctree), horiz=T)
dev.off()

# linear regression for location_name with death perc
cause_death_model=glm(perc~cause_name, data=perc_dt)
coeff_table=summary(cause_death_model)$coeff[-1,]
rownames(coeff_table)=gsub("cause_name", '', rownames(coeff_table))
sig_cause=coeff_table[coeff_table[,4]<0.05,]
sig_cause=sig_cause[order(sig_cause[,4]),]
top10_causes=sig_cause[1:10,]   # Most common top 10 causes in each country, not top10 for the world

# select out significant cause of death
top10_causes_dt=subset(perc_dt, cause_name %in% rownames(top10_causes))



# Heatmap for top10 cause in different world region
top10_cause_by_region=aggregate(top10_causes_dt, 
            by=list(top10_causes_dt$region,
                    top10_causes_dt$cause_name), FUN=mean)[,c('Group.1', 'Group.2', 'perc')]
colnames(top10_cause_by_region)=c("region",'cause_name', 'perc')
top10_cause_by_region$cause_name<- factor(top10_cause_by_region$cause_name, # order the cause in plot based on top 10 order 
                      levels=rownames(top10_causes))
  
top10_region_heatmap=ggplot(top10_cause_by_region, aes(cause_name, region)) +
   geom_tile(aes(fill = perc), colour = "white") +    # tile setup
   scale_fill_gradient(low = "white", high="blue")+   # tile coloring
   theme(axis.text.x=element_text(angle=60, hjust=1)) # rotate x-label to vertical

# Save heatmap into a file in png format
ggsave(file="results/top10_region_heatmap.png", 
       width=4, height=5, units='in', # define dimension of the png 
       plot=top10_region_heatmap)

# Save important variables for future analysis
save(cause_death_model,
     top10_causes,
     top10_cause_by_region, 
     file='results/significant_cause_of_death.Rdata')
