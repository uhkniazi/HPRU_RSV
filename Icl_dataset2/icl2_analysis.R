# File: icl2_analysis.R
# Desc: Data set from duke uni, affy
# Auth: umar.niazi@kcl.ac.uk
# Date: 07/10/2016

gcwd = getwd()
setwd('Icl_dataset2/')
source('Header.R')

#### data loading
##### connect to mysql database to get samples
library(RMySQL)
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
# another way to get the query, preferred
dfSample = dbGetQuery(db, "select * from MetaFile where idData=8")
# close connection after getting data
dbDisconnect(db)

n = paste0(dfSample$location[dfSample$id == 20], dfSample$name[dfSample$id == 20])
load(n)
oExp.lumi = lumi.n.q
rm(lumi.n.q)

## data formatting for DE analysis
i = which(oExp.lumi$Symptoms == 'Y')
oExp.lumi = oExp.lumi[,i]
# keep only days 0 to 7
i = which(oExp.lumi$day %in% c(0, 1, 3, 7, 10))
oExp.lumi = oExp.lumi[,i]
mDat = exprs(oExp.lumi)
# sample annotation
dfSamples = pData(oExp.lumi)
dfSamples = droplevels.data.frame(dfSamples)

# sanity check
identical(as.character(dfSamples$samplesNum), colnames(mDat))

xtabs( ~ PID + day, data=dfSamples)

## create grouping factors by days
fSamples = factor(oExp.lumi$day)
table(fSamples)
levels(fSamples)


# select data matrix and remove genes with low detection call
ivDetection = detectionCall(oExp.lumi)
mDat = mDat[ivDetection > 0, ]
####### end grouping

################################### main analysis section
# map the nuID to human symbols
cvSym = getSYMBOL(rownames(mDat), 'lumiHumanAll.db')
# number of genes not annotated
table(is.na(cvSym))
# remove unannotated genes
mDat = mDat[!is.na(cvSym),]

### lme library and function
library(lme4)
library(car)

f_get.lme.pvalue = function(x, fac, ran){
  return(Anova(lmer(x ~ fac + (1|ran), REML = F))[[3]])
}

ivSigGenes = apply(mDat, 1, f_get.lme.pvalue, fSamples, dfSamples$PID)
ivSigGenes.adj = p.adjust(ivSigGenes, 'BH')

# select the genes with significant model deviance p-values
i = which(ivSigGenes.adj < 0.01)

mDat = mDat[i,]

# load library to assign p-values to coefficients
library(lmerTest)
# refit the models on these significant genes 
lFit = lapply(1:nrow(mDat), function(x){
  lmer(mDat[x,] ~ fSamples + (1|dfSamples$PID), REML=F)
})

names(lFit) = rownames(mDat)

f_get.coef.pvalue = function(fit){
  s = summary(fit)
  return(s$coefficients[,'Pr(>|t|)'])
}

mFit.pval = sapply(lFit, f_get.coef.pvalue)
rownames(mFit.pval) = levels(fSamples)

dfGenes = select(lumiHumanAll.db, keys = rownames(mDat), columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# sanity check
nrow(mDat) == nrow(dfGenes)
rownames(dfGenes) = dfGenes$PROBEID

# get the list of genes for each comparison i.e. each coefficient compared to base line
lSigGenes.adj = vector('list', length = length(levels(fSamples))-1)
names(lSigGenes.adj) = levels(fSamples)[2:length(levels(fSamples))]

for (i in 2:length(levels(fSamples))){
  p.val = mFit.pval[i,]
  lSigGenes.adj[[i-1]] = names(p.val)[p.val < 0.05]
}

sapply(lSigGenes.adj, length)

# grouping of genes based on expression trends
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
cvCommonGenes = NULL
for (i in seq_along(n)) {
  cvCommonGenes = append(cvCommonGenes, lSigGenes.adj[[names(n[i])]])
}
cvCommonGenes = unique(cvCommonGenes)
mCommonGenes = sapply(seq_along(n), function(x) cvCommonGenes %in% lSigGenes.adj[[names(n[x])]])
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(n)

#### analysis by grouping genes
# create groups in the data based on 2^n-1 combinations
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)
plot(hc)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.7)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mCommonGenes.grp = cbind(mCommonGenes.grp, cp)

rn = apply(mCommonGenes, 1, function(x) any(x[c(1,2)] == c(T, T)))
rn = names(rn[rn])
df.rn = select(lumiHumanAll.db, keys = rn, columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# write csv to look at gene list
write.csv(df.rn[,-1], file=paste('Temp/', 'early_response_commongenes', '.csv', sep=''))

# mid response phase days 3 to 5
rn = apply(mCommonGenes, 1, function(x) any(x[c(2, 3)] == c(T, T)))
rn = names(rn[rn])

df.rn = select(lumiHumanAll.db, keys = rn, columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# write csv to look at gene list
write.csv(df.rn[,-1], file=paste('Temp/', 'mid_response_commongenes', '.csv', sep=''))

# late response genes 5 to 7
rn = apply(mCommonGenes, 1, function(x) any(x[c(3, 4)] == c(T, T)))
rn = names(rn[rn])
df.rn = select(lumiHumanAll.db, keys = rn, columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# write csv to look at gene list
write.csv(df.rn[,-1], file=paste('Temp/', 'late_response_commongenes', '.csv', sep=''))

## all overexpressed genes if interested in
i = 1:nrow(mCommonGenes)

m1 = as.matrix(mCommonGenes[i,])
m1 = mDat[rownames(m1),]

fGroups = fSamples
colnames(m1) = fGroups
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]

## download the cgraph library
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/master/CGraphClust.R'
download(url, 'CGraphClust.R')

# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')
library('NMF')
# ignore step if stabalization not required
m1 = t(apply(m1, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(m1) = fGroups

# scale across rows
rownames(m1) = dfGenes[rownames(m1), 'SYMBOL']
mCounts = t(m1)
mCounts = scale(mCounts)
mCounts = t(mCounts)
# threshhold the values
mCounts[mCounts < -3] = -3
mCounts[mCounts > 3] = 3

# draw the heatmap  color='-RdBu:50'
aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)


########### pathway analysis using CGraph library
# uniprot annotation for data
cvGenes = rownames(mCommonGenes)
dfGenes = select(lumiHumanAll.db, keys = cvGenes, columns = c('ENTREZID', 'SYMBOL', 'GENENAME', 'UNIPROT'), keytype = 'PROBEID')
dfGenes = na.omit(dfGenes)

# get reactome data
url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
dir.create('Data_external', showWarnings = F)
csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
# map reactome pathways
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
x = gsub('\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfGenes$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfGenes$UNIPROT)
dfReactome.sub$ENTREZID = dfGenes$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)

# get expression data
mCounts = mDat[unique(dfGenes$PROBEID),]
## optional adjusting the data for repeated measurements
mCounts = sapply(rownames(mCounts), function(x){
  return(fitted(lFit[[x]]))
})
mCounts = t(mCounts)
# map the probe names to enterez ids
i = match(rownames(mCounts), dfGenes$PROBEID)
rownames(mCounts) = dfGenes$ENTREZID[i]
fGroups = fSamples
names(fGroups) = dfSamples$PID
colnames(mCounts) = fGroups
# reorder on grouping factor
mCounts = mCounts[,order(fGroups)]
fGroups = fGroups[order(fGroups)]
mCounts = t(mCounts)

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
levels(fGroups)

# create a correlation matrix to decide cor cutoff
mCor = cor(mCounts)

# check distribution 
hist(mCor, prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.5, bSuppressPlots = F)

## general graph structure
set.seed(1)
plot.final.graph(oGr)
ecount(getFinalGraph(oGr))
vcount(getFinalGraph(oGr))
## community structure
# plot the main communities and centrality graphs
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 40)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
set.seed(1)
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 20)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')

# get community sizes
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
iSizes = sort(table(dfCluster$cluster))
# remove communities smaller than 5 members or choose a size of your liking
i = which(iSizes <= 5)
if (length(i) > 0) {
  cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
  iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
  oGr = CGraphClust.recalibrate(oGr, iVertKeep)
}

## centrality diagnostics
## centrality parameters should not be correlated significantly and the location of the central
## genes can be visualized
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
# plot the genes or vertex sizes by fold change
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)

## the diagnostic plots show the distribution of the centrality parameters
# these diagnostics plots should be looked at in combination with the centrality graphs
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)

## top vertices based on centrality scores
## get a table of top vertices 
dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.90)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)

####### NOTE: This section of code is very slow, use only if you need data from genbank
# loop and get the data from genbank
n = rep(NA, length=nrow(dfTopGenes.cent))
names(n) = as.character(dfTopGenes.cent$VertexID)
for (i in seq_along(n)){
  n[i] = f_csGetGeneSummaryFromGenbank(names(n)[i])
  # this wait time is required to stop sending queries to ncbi server very quickly
  Sys.sleep(time = 3)
}
cvSum.2 = as.character(dfTopGenes.cent$VertexID)
dfTopGenes.cent$Summary = n[cvSum.2]
####### Section ends

write.csv(dfTopGenes.cent, file='Temp/Top_Centrality_Genes_icl2.csv')

## if we want to look at the expression profiles of the top genes
# plot a heatmap of these top genes
library(NMF)
m1 = mCounts[,as.character(dfTopGenes.cent$VertexID)]
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)

m1 = mCounts[,as.character(dfTopGenes.cent$VertexID)]
m1 = apply(m1, 2, f_ivStabilizeData, fGroups)
rownames(m1) = fGroups
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)

## in addition to heatmaps the graphs can be plotted
# plot a graph of these top genes
# plot for each contrast i.e. base line vs other level
lev = levels(fGroups)[-1]
m = mCounts
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = induced_subgraph(getFinalGraph(oGr), vids = as.character(dfTopGenes.cent$VertexID))
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=50)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.2, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}


### Looking at the largest clique can be informative in the graph
# plot the graph with location of the clique highlighted
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique at each grouping contrast
lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=50)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, layout=layout_with_fr, main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}



## instead of looking at individual genes we can look at clusters
## we can look at the problem from the other direction and look at clusters instead of genes
# some sample plots
# mean expression of groups in every cluster
par(p.old)
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster', cex.axis=0.7)
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)#, p.cut=0.05 )
# principal component plots
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)#, p.cut=0.05)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)#, p.cut=0.05)
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)#, p.cut=0.05)
# plot variance of cluster
m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters

csClust = rownames(m)
length(csClust)
pdf('Temp/cluster_variance_duke.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(m[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())
#boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust))

# print the labels of the clusters from reactome table
i = which(dfReactome.sub$V2 %in% csClust)
dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
rownames(dfCluster.name) = NULL
dfCluster.name

#### plot a graph of clusters 
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
data.frame(sort(table(dfCluster$cluster)))
#csClust = rownames(m$clusters)
csClust = as.character(unique(dfCluster$cluster))

mMarginal = getClusterMarginal(oGr, t(mCounts))

library(lattice)
dfData = data.frame(t(mMarginal))
dfData$days = fGroups
dfData$id = factor(names(fGroups))

xyplot(X1280218 ~ days | id, data=dfData, type='o')
n = colnames(dfData)[1:15]

sapply(n, function(x) {p =xyplot(dfData[,x] ~ days | id, data=dfData, type='o', ylab=x)
print(p)})

sapply(n, function(x) {p =bwplot(dfData[,x] ~ days, data=dfData, type='o', ylab=x)
print(p)})

dfStack = stack(dfData)
dfStack$days = dfData$days

xyplot(values ~ days | ind, data=dfStack, type='o')

## xyplots with population average
mMarginal = getClusterMarginal(oGr, t(mCounts))
mMarginal = apply(mMarginal, 1, function(x) tapply(x, fGroups, mean))
mMarginal = scale(mMarginal)
dfData = data.frame(mMarginal)
cn = colnames(mMarginal)
rownames(dfCluster.name) = dfCluster.name$V2
cn = c('Interleukins Signalling', 'Adaptive Immunity', 'Membrane Traffic', 'Neutrophil Degranulation', 'Interferon Signalling',
       'GPCR Signalling', 'Lipids metabolism', 'Protein metabolism')
colnames(dfData) = cn
#colnames(dfData) = gsub('X', '', colnames(dfData))
dfData$days = factor(rownames(dfData), levels=c(0, 1, 3, 7, 10))

n = colnames(dfData)[1:8]

sapply(n, function(x) {p =xyplot(dfData[,x] ~ days , data=dfData, type='o', ylab=x)
print(p)})

dfStack = stack(dfData)
dfStack$days = dfData$days

xyplot(values ~ days | ind, data=dfStack, type='o')





