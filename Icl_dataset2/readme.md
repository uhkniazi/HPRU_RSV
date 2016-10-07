Import and analysis of the ICL RSV Challenge dataset from ICL.

# Database entries  
Project.ID = 3  
Data.id = 8 (ICL second dataset)  
[Evernote link](https://www.evernote.com/shard/s288/nl/38698211/9dd55b6b-87ed-44eb-827c-96e897e8fb0a?title=00%20HPRU%20RSV "Private")

# Scripts
**In no particular order**  

1. **Header.R**  
  * General header file for setting variables and functions.  
2. **icl2_import.R**  
  * Import the illumina raw data, normalise the data, create appropriate database entries and save the expression set object.  
3. **icl2_analysis.R**  
  * Load the lumi object and relevant covariates, fit a linear mixed model to each gene, and select genes that show significant DE. Cluster coexpressed genes using CGraphClust library.  
  

