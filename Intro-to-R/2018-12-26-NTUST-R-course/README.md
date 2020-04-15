# Introduction to R

This repository provides simple introdcution to R

## Basic Data Analysis

This section performs clustering analysis on a simple cause of death data set to identify top10 cause of death by country and summarize by region.  

1. Data Set: IHME_GBD_2017_PreProcessed_DataSets.Rdata
    - Downloaded from global data exchange (<http://ghdx.healthdata.org/>)
    - Contain all cause of death from 195 countries for 2017
    - Data is represented by the percentage of cause of death by country
    - Raw data is stored in data/IHME-GBD/
    - Raw data was processed by data_preprocessing.R 
    - Processed data was stored as IHME_GBD_2017_PreProcessed_DataSets.Rdata object

2. Script: cluster_code.R
    - The script loads reshape2, tidyverse libraries and Data Set .Rdata object
    - Create a 'results' subdirectory if it does exists
    - Performs PCA, store image in results
    - Perform Hierarchical Clustering, store image in results
    - Perform simple general linear model to select top 10 death causes
    - Summarize the top10 causes by geographical regions and generate a heatmap
    - Save the heatmap




| | B | C | D
|---|---|---|---
A|1|2|3
B|1|2|3
C|1|2|3
