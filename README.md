# HPRU_RSV
RSV challenge data analysis

# HJ_Habibi_200213_Apr13_No.R
Dataset 1, /home/uniazi/Data/R/HPRU_RSV/Data_external/Original/RSV_Challenge_for_Umar/RSV_Expt 1_Training Set_parameters.csv  
The following steps are done for analysis:  
1- Data loading, formatting, normalization using lumi.  
2- Removing outliers using PCA.  
3- Removing genes with 0 detection calls and no annotation.  
4- Selecting factors, creating design matrix, limma to call DE.  
5- Selecting Differentially expressed genes in each factor level compared to base line.  
6- Producing volcano plots, heatmaps, csv files and a graph. GO term enrichment analysis done using innate db.  

# HJ_Habibi_200213_Apr13_No_Grps_1_2_3.R
Same dataset as HJ_Habibi_200213_Apr13_No.R, however the grouping factor was Study Group 1, 2 and 3.  
The analysis was similar to the previous script. However as there were more DE genes in this case, a more detailed CGraph based 
clustering was performed, with various summary plots, particularly the top genes based on centrality scores.

# HJ_Habibi_200213_Apr13_No_Symptoms.R
Same dataset as HJ_Habibi_200213_Apr13_No.R, however the grouping factor was Symptoms based on viral load. 



