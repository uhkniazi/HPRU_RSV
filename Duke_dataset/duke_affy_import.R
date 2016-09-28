# File: duke_affy_import.R
# Desc: Data set from duke uni, affy
# Auth: umar.niazi@kcl.ac.uk
# Date: 28/09/2016


gcwd = getwd()
setwd('Duke_dataset/')
source('Header.R')


#### data loading
cvDir = getwd()
setwd('Data_external/Duke_data')
setwd('2013_rsv_cel/')
# set working directory to the directory with all affymetrix cel files
oData = ReadAffy()
# normalise the data
x.affy = rma(oData)

# check the probe annotation file
annotation(x.affy)
library(hgu133a2.db)

setwd(cvDir)
## load the sample annotation data set
dfSamples = read.csv(file.choose(), header=T, stringsAsFactors=F)

rn = rownames(pData(x.affy))
# order the sample names according to cel files
table(rn %in% dfSamples$cel_name)
i = match(rn, dfSamples$cel_name)
dfSamples = dfSamples[i,]

# sanity check
identical(rn, dfSamples$cel_name)

pData(x.affy) = dfSamples

# load symptoms data to set the covariates/predictors for the samples
df = read.csv(file.choose(), header=T, stringsAsFactors=F)
# the infected samples that showed symptoms only
table(dfSamples$subject_id %in% df$PID)

f = dfSamples$subject_id %in% df$PID
# remove this missing sample from array object
x.affy = x.affy[,f]
dfSample = pData(x.affy)
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
cn = dbListFields(db, 'Sample')[-1]
df = data.frame(3, 5, dfSample$subject_id, dfSample$cel_name, 'RSV Dataset Affymetrix Sample from Duke University for infected symptomatic at particular days',
                dfSample$timepoint, dfSample$sample_id, 'Infected Symptomatic')
colnames(df) = cn
dbWriteTable(db, name = 'Sample', value=df, append=T, row.names=F)

# save the affymetrix object
dbListTables(db)
cn = dbListFields(db, 'MetaFile')[-1]
n = make.names('oAffymetrixExpressionSet for RSV Duke Infected Symptomatic rds')
df = data.frame(5, n, '.rds', '~/Data/MetaData/', 'RSV Dataset Affymetrix Sample from Duke University for infected symptomatic at particular days')
colnames(df) = cn
n2 = paste0('Objects/', n)
save(x.affy, file=n2)
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)



