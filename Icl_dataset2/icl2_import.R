# File: icl2_import.R
# Desc: Data set from icl, illumina
# Auth: umar.niazi@kcl.ac.uk
# Date: 07/10/2016

gcwd = getwd()
setwd('Icl_dataset2/')
source('Header.R')

#### data loading
## load the annotation data set
dfSamples = read.csv('Data_external/annotation.csv', header=T)
# load the second annotation of outcomes
dfOutcomes = read.csv('Data_external/outcome.csv', header=T)
# match the two data frames
i = match(dfSamples$subjectID, dfOutcomes$PID)

df = data.frame(samplesNum=dfSamples$samplesNum, PID=dfOutcomes$PID[i], day=dfSamples$day, Outcome=dfOutcomes$Outcome[i],
                Symptoms=dfOutcomes$Symptoms[i])
# remove one missing row from the end
dfSamples = df[1:(nrow(df)-1),]
dfSamples = droplevels.data.frame(dfSamples)
str(dfSamples)
rownames(dfSamples) = dfSamples$samplesNum


# load the raw data
csFile.1 = 'Data_external/raw_data.txt'

# create lumi object
x.lumi = lumiR.batch(c(csFile.1), lib.mapping = 'lumiHumanIDMapping')

# add the sample annotation
# check if sample ids same
t1 = colnames(x.lumi)
t2 = rownames(dfSamples)
table(t1 %in% t2)
# remove missing samples
f = t1 %in% t2
x.lumi = x.lumi[,f]

# sanity check
t1 = colnames(x.lumi)
t2 = rownames(dfSamples)
table(t1 %in% t2)

df = pData(x.lumi)
# sort in same order
dfSamples = dfSamples[rownames(df),]
dfSamples$sampleID = df$sampleID
dfSamples = droplevels.data.frame(dfSamples)
pData(x.lumi) = dfSamples
rm(dfSamples)

#### initial quality checks using PCA
m = exprs(x.lumi)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(x.lumi$day)
fSamples = as.factor(x.lumi$Outcome)

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3, not normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3, not normalized')
par(p.old)


#################### variance stabalization and normalization
lumi.t = lumiT(x.lumi)
lumi.n = lumiN(lumi.t, method = 'rsn')

################## QC checks after normalization
lumi.n.q = lumiQ(lumi.n)
plot(lumi.n.q)

m = exprs(lumi.n.q)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(lumi.n.q$day)
fSamples = as.factor(lumi.n.q$Outcome)

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], ncol=2)
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')
par(p.old)

# remove outlier
l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)
# i = which(l$cluster.label$c1 %in% c(1))
# lumi.n.q = lumi.n.q[,-i]
# 
# m = exprs(lumi.n.q)
# # pca on samples i.e. covariance matrix of m
# pr.out = prcomp(t(m), scale=T)
# ## choose appropriate factor
# fSamples = as.factor(lumi.n.q$day)
# fSamples = as.factor(lumi.n.q$Outcome)
# fSamples = lumi.n.q$PID
# 
# col.p = rainbow(length(unique(fSamples)))
# col = col.p[as.numeric(fSamples)]
# # plot the pca components
# par(mfrow=c(2,2))
# plot.new()
# legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], ncol=2)
# plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
#      main='PCA comp 1 and 2')
# plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
#      main='PCA comp 1 and 3')
# plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
#      main='PCA comp 2 and 3')
# par(p.old)

## create the database entries
dfSample = pData(lumi.n.q)
library('RMySQL')
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
cn = dbListFields(db, 'Sample')[-1]
df = data.frame(3, 8, dfSample$PID, dfSample$samplesNum, 
                'RSV Dataset Illumina microarray Sample from Imperial College London for infected symptomatic and asymptomatic at particular days, covariates: timepoint, infection status, symptoms',
                dfSample$day, dfSample$Outcome, dfSample$Symptoms)
colnames(df) = cn
dbWriteTable(db, name = 'Sample', value=df, append=T, row.names=F)

# save the lumi object
dbListTables(db)
cn = dbListFields(db, 'MetaFile')[-1]
n = make.names('oLumiExpressionSet for RSV ICL Sep 2013 rds')
df = data.frame(8, n, '.rds', '~/Data/MetaData/', 'RSV Dataset Illumina microarray Sample from Imperial College London for infected symptomatic and asymptomatic at particular days, covariates: timepoint, infection status, symptoms')
colnames(df) = cn
n2 = paste0('Objects/', n)
save(lumi.n.q, file=n2)
dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
dbDisconnect(db)


