library(reshape2)
library(tidyverse)
library(readxl)
#
# read in raw data table in 'melt' form and 
# convert the table into a matrix of percentage
# 

raw_dt=read.csv('data/IHME-GBD_2017_DATA/IHME_GBD_2017_DATA.csv')
location_code=read_xlsx('data/IHME-GBD_2017_DATA/IHME_GBD_2017_GBD_LOCATIONS_HIERARCHY.XLSX')
raw_location=merge(raw_dt, location_code, by='location_id')
raw_location=merge(raw_location, location_code,
                   by.x='parent_id',
                   by.y='location_id')

perc_dt=subset(raw_location, metric_id==2, select=c('measure_name', 
                                              'location_name.x',
                                              'location_name',
                                              'sex_name',
                                              'cause_name',
                                              'metric_name',
                                              'year',
                                              'val'))
colnames(perc_dt)[2]='country'
colnames(perc_dt)[3]='region'
colnames(perc_dt)[8]='perc'

matrix_dt=dcast(perc_dt, cause_name~country, value.var="perc")
rownames(matrix_dt)=matrix_dt[,1]
matrix_dt=t(matrix_dt[,-1])
regions=lapply(unique(perc_dt$region), function(X, region_data){
  unique(subset(region_data, region==X, country))
  }, region_data=perc_dt)
names(regions)=unique(perc_dt$region)
  

save(matrix_dt, perc_dt, raw_dt, location_code, regions, 
     file='data/IHME_GBD_2017_PreProcessed_DataSets.Rdata')

# END