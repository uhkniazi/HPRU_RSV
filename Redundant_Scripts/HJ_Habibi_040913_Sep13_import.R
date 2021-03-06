# File: HJ_Habibi_040913_Sep13_import.R
# Desc: Data set from RSV challenge - part 2
# Auth: u.niazi@imperial.ac.uk
# Date: 03/02/16

source('Header.R')

#### data loading
## load the first data set
dfSamples = read.csv(file.choose(), header=T)
rownames(dfSamples) = dfSamples$samplesNum
csFile = 'Data_external/Original/RSV_Challenge_for_Umar/RSV_challange_handover(Artem)/original annotation and data files/HJ_Habibi_040913_NouRNA_noNorn_bkgsubt__Sample_Probe_Profile.txt'

# create lumi object
x.lumi = lumiR.batch(csFile, lib.mapping = 'lumiHumanIDMapping')

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

#### initial quality checks using PCA
m = exprs(x.lumi)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(x.lumi$day)
fSamples = as.factor(x.lumi$group)
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

#################### variance stabalization and normalization
lumi.t = lumiT(x.lumi)
lumi.n = lumiN(lumi.t, method = 'rsn')

################## QC checks after normalization
lumi.n.q = lumiQ(lumi.n)

m = exprs(lumi.n.q)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(lumi.n.q$day)
fSamples = as.factor(as.character(lumi.n.q$ethnicity))
fSamples = as.factor(lumi.n.q$group)

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
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

# # remove the outliers 
# # remove the outlier groups from the data
# # these can be seen on the pc2 and pc3 plots
# m = pr.out$x[,1:3]
# m = data.frame(m, fSamples)
# i = which(m$PC2 > 0 & m$PC1 < -300)
# i = unique(c(i, which(m$PC3 < -50)))
# c = col
# c[i] = 'black'
# par(mfrow=c(2,2))
# plot.new()
# legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
# plot(pr.out$x[,1:2], col=c, pch=19, xlab='Z1', ylab='Z2',
#      main='PCA comp 1 and 2')
# plot(pr.out$x[,c(1,3)], col=c, pch=19, xlab='Z1', ylab='Z3',
#      main='PCA comp 1 and 3')
# plot(pr.out$x[,c(2,3)], col=c, pch=19, xlab='Z2', ylab='Z3',
#      main='PCA comp 2 and 3')
# par(p.old)
# # identify these groups
# pData(lumi.n.q)[i,]

# remove the outliers
#oExp.lumi = lumi.n.q[,-i]
oExp.lumi = lumi.n.q

# # plot another PCA
# m = exprs(oExp.lumi)
# # pca on samples i.e. covariance matrix of m
# pr.out = prcomp(t(m), scale=T)
# ## choose appropriate factor
# fSamples = as.factor(oExp.lumi$Day)
# fSamples = as.factor(oExp.lumi$Sample.group.3)
# fSamples = as.factor(oExp.lumi$Study.Group)
# 
# col.p = rainbow(length(unique(fSamples)))
# col = col.p[as.numeric(fSamples)]
# # plot the pca components
# par(mfrow=c(2,2))
# plot.new()
# legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
# plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
#      main='PCA comp 1 and 2')
# plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
#      main='PCA comp 1 and 3')
# plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
#      main='PCA comp 2 and 3')
# par(p.old)
# f_Plot3DPCA(pr.out$x[,1:3], col, pch=19, xlab='Z1', ylab='Z2', zlab='Z3', 
#             main='Plot of first 3 components')
# par(p.old)

# save the normalized object
save(oExp.lumi, file='Objects/lumi.n.experiment2_HJ_Habibi_040913_Sep13.rds')

