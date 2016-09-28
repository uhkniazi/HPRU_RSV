# File: duke_affy_import.R
# Desc: Data set from duke uni, affy
# Auth: umar.niazi@kcl.ac.uk
# Date: 28/09/2016

source('Header.R')
library(affy)


#### data loading
cvDir = getwd()
setwd('Data_external/2013_RSV/2013_rsv_cel/')
oData = ReadAffy()
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

# load symptoms data
df = read.csv(file.choose(), header=T, stringsAsFactors=F)
i = which(df$sample_id == '')
df = df[-i,]
table(dfSamples$sample_id %in% df$sample_id)
i = dfSamples$sample_id %in% df$sample_id
dfSamples = dfSamples[i,]
# remove this missing sample from array object
x.affy = x.affy[,i]

# some samples do not match the ids in the symptoms table
table(df$sample_id %in% dfSamples$sample_id)
# order symptoms table accordingly
i = match(dfSamples$sample_id, df$sample_id)
df = df[i,]

# load serology outcome data
df = read.csv('Data_external/Serology data/outcome.csv', header=T, stringsAsFactors=F)
# choose patients whose ids match with loaded samples
i = df$PID %in% dfSamples$subject_id
df = df[i,]
# order the data frame accordingly
table(dfSamples$subject_id %in% df$PID)
# remove non-matching cel and patient ids
i = (dfSamples$subject_id %in% df$PID)
dfSamples = dfSamples[i,]
# remove this missing sample from array object
x.affy = x.affy[,i]
# order df according to patient ids
i = match(dfSamples$subject_id, df$PID)
df = df[i,]
dfSamples$outcome = df$Outcome
dfSamples$symptoms = df$Symptoms

## add this information to the expression set object
rownames(dfSamples) = dfSamples$cel_name

# add the sample annotation
# check if sample ids same
t1 = colnames(x.affy)
t2 = rownames(dfSamples)
table(t1 %in% t2)
# remove missing samples
f = t1 %in% t2
x.affy = x.affy[,f]

# sanity check
t1 = colnames(x.affy)
t2 = rownames(dfSamples)
table(t1 %in% t2)

pData(x.affy) = dfSamples
rm(dfSamples)
#### initial quality checks using PCA
m = exprs(x.affy)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(x.affy$outcome)
fSamples = as.factor(x.affy$timepoint)

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], ncol = 2)
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3, not normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3, not normalized')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components - not normalized')
par(p.old)

save(x.affy, file='Objects/duke_normalised_inf_uninf.rds')

## remove uninfected
i = x.affy$outcome
i = which(i == 'Infected')
x.affy.all = x.affy
x.affy = x.affy[,i]

## pca plots to check for outliers
m = exprs(x.affy)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(x.affy$timepoint)

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], ncol = 2)
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3, not normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3, not normalized')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components - not normalized')
par(p.old)

l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)

## add the combat step 
fBatch = factor(x.affy$subject_id)

library(sva)
modcombat = model.matrix(~1, data=pData(x.affy))
oCexp = ComBat(exprs(x.affy), batch = fBatch, mod=modcombat)
m = oCexp
pr.out = prcomp(t(m), scale=T)
fSamples = as.factor(x.affy$timepoint)
fSamples = fBatch
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], ncol = 2)
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2, not normalized')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3, not normalized')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3, not normalized')
par(p.old)

exprs(x.affy) = oCexp

# remove outliers
l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)
i = which(l$cluster.label$c1 == '1')
i2 = which(l$cluster.label$c1 == '9')
i = unique(c(i, i2))
# remove the outliers 
c = col
c[i] = 'black'
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
plot(pr.out$x[,1:2], col=c, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=c, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=c, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')
par(p.old)
# identify these groups
pData(x.affy)[i,]

# remove the outliers
x.affy = x.affy[,-i]

# plot another PCA
m = exprs(x.affy)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(x.affy$timepoint)

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))], ncol = 2)
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')
par(p.old)
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components')
par(p.old)

l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)

# save the normalized object
save(x.affy, file='Objects/duke_normalised_inf_combat_individuals.rds')

