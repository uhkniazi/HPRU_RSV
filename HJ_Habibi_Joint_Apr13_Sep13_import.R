# File: HJ_Habibi_Joint_Apr13_Sep13_import.R
# Desc: Data set from RSV challenge - joint analysis of both datasets
# Auth: u.niazi@imperial.ac.uk
# Date: 03/02/16

source('Header.R')

#### data loading
## load the annotation data set
dfSamples = read.csv('Data_external/Original/RSV_Challenge_for_Umar/RSV_challenge_annotation_MH.csv', header=T)
rownames(dfSamples) = dfSamples$samplesNum
head(dfSamples)

# load the 2 datasets
csFile.1 = 'Data_external/Original/RSV_Challenge_for_Umar/RSV_challange_handover(Artem)/original annotation and data files/HJ_Habibi_040913_NouRNA_noNorn_bkgsubt__Sample_Probe_Profile.txt'
csFile.2 = 'Data_external/Original/RSV_Challenge_for_Umar/RSV_challange_handover(Artem)/original annotation and data files/HJ_Habibi_200213_Apr13_No uRNA_Bkg sub_No Norm_Sample_Probe_Profile.txt'

# create lumi object
x.lumi = lumiR.batch(c(csFile.1, csFile.2), lib.mapping = 'lumiHumanIDMapping')

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

pData(x.lumi) = dfSamples
rm(dfSamples)

# create a batch variable
m = expand.grid(0:1, 2:9)
fBatch = paste0('MH0', m[,1], m[,2])
fBatch = append(fBatch, 'MH022')
i = which(x.lumi$subjectID %in% fBatch)
# create two batches 1 and 2 representing two experiments
fBatch = rep(2, length.out = ncol(x.lumi))
fBatch[i] = 1
# sanity check
table(fBatch)
x.lumi$fBatch = factor(fBatch)
#### initial quality checks using PCA
m = exprs(x.lumi)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(x.lumi$day)
fSamples = as.factor(x.lumi$group)
fSamples = as.factor(x.lumi$fBatch)
#fSamples = as.factor(x.lumi$Study.Group)

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
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components - not normalized')
par(p.old)

l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)

i = which(l$cluster.label$c1 %in% c(1:4))
fBatch.2 = rep(2, length.out = ncol(x.lumi))
fBatch.2[i] = 1
x.lumi$fBatch.2 = factor(fBatch.2)
fSamples = x.lumi$fBatch.2
# check again

#################### variance stabalization and normalization
lumi.t = lumiT(x.lumi)
## add the combat step 
library(sva)
modcombat = model.matrix(~ 1, data=pData(lumi.t))
oCexp = ComBat(exprs(lumi.t), batch = lumi.t$fBatch.2, mod=modcombat)
exprs(lumi.t) = oCexp

lumi.n = lumiN(lumi.t, method = 'rsn')

################## QC checks after normalization
lumi.n.q = lumiQ(lumi.n)

m = exprs(lumi.n.q)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(lumi.n.q$day)
# fSamples = as.factor(as.character(lumi.n.q$ethnicity))
# fSamples = as.factor(lumi.n.q$group)
# fSamples = as.factor(lumi.n.q$fBatch)
fSamples = as.factor(lumi.n.q$fBatch.2)
fSamples = as.factor(lumi.n.q$subjectID)

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
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components')
par(p.old)

# remove outlier
l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)
i = which(l$cluster.label$c1 %in% c(1:3))
lumi.n.q = lumi.n.q[,-i]

## add the combat step 
library(sva)
modcombat = model.matrix(~ 1, data=pData(lumi.n.q))
oCexp = ComBat(exprs(lumi.n.q), batch = factor(lumi.n.q$subjectID), mod=modcombat)
exprs(lumi.n.q) = oCexp

## check after combat and normalization
m = exprs(lumi.n.q)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(lumi.n.q$day)
# fSamples = as.factor(as.character(lumi.n.q$ethnicity))
# fSamples = as.factor(lumi.n.q$group)
# fSamples = as.factor(lumi.n.q$fBatch.3)

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
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components')
par(p.old)

l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)

i = which(l$cluster.label$c1 %in% 1:2)
# remove the outliers 
# remove the outlier groups from the data
# these can be seen on the pc2 and pc3 plots
# m = pr.out$x[,1:3]
# m = data.frame(m, fSamples)
# i = which(m$PC2 > 0 & m$PC1 < -300)
# i = unique(c(i, which(m$PC3 < -50)))
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
#pData(lumi.n.q)[i,]

# remove the outliers
oExp.lumi = lumi.n.q[,-i]

# plot another PCA
m = exprs(oExp.lumi)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(oExp.lumi$group)

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
f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
            main='Plot of first 3 components')
par(p.old)

l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)

# save the normalized object
save(oExp.lumi, file='Objects/lumi.n.Combined_Apr13_Sep13_adjusted.rds')

