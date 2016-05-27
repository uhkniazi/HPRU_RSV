# HPRU_RSV
RSV challenge data analysis

# HJ_Habibi_040913_Sep13_Analysis1_2.R
Sep 2013 dataset

# HJ_Habibi_Joint_Apr13_Sep13_import.R
Joint analysis of both april and september datasets, merging using Combat.
## HJ_Habibi_Joint_Apr13_Sep13_Analysis1.R
Main analysis after creating joint lumi object

# viral_load_clustering.R, antibody_clustering.R
Scripts to cluster the data based on viral load or antibody responses, in order to define groupings for the microarray analysis.  
https://www.evernote.com/shard/s288/nl/38698211/e9c026fa-1e40-4389-814f-b9e5a1bc50e8?title=New%20Groupings%20with%20PCR%20Results

# antibody_clustering_inf_uninf.R
extension of the antibody_clustering.R script with both infected and uninfected groups, and creating AB1, AB2, AB3 groups.  
https://www.evernote.com/shard/s288/nl/38698211/16a94625-014e-46d2-96dc-eef8872ed29c?title=Antibody%20clustering%2C%20infected%20and%20uninfected  



# Redundant_Scripts
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


# aggie_pcr.R
Analysis of some qpcr ct values, and generating missing data using a posterior predictive distribution and non-informative prior.

# HJ_Habibi_200213_Apr13_Analysis10.R and Analysis11
Both the analyses are similar, but we remove asymptomatic infected in analysis11 and do a joint analysis.  
see https://www.evernote.com/shard/s288/nl/38698211/b308257b-2786-4e94-8b7e-485ac12272e2 and end of report for analysis 11.

# HJ_Habibi_040913_Sep13_Analysis1_2.R
Analysis 1 for the second data set. Details here: https://www.evernote.com/shard/s288/nl/38698211/43a5359b-fce7-4dad-8b92-36a75e4b4623

# HJ_Habibi_Joint_Apr13_Sep13_import.R
import and merging of the 2 data sets using combat

# HJ_Habibi_Joint_Apr13_Sep13_Analysis1.R
analysis of the joint data set, with graph analysis and some variable selection for class prediction and roc curves.

# viral_load_clustering.R, antibody_clustering.R
various unsupervised clustering methods to identify groups in the data.  
