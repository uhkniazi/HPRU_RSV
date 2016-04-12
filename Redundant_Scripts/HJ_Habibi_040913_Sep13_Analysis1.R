# File: HJ_Habibi_040913_Sep13_Analysis1.R
# Desc: Data set from RSV challenge experiment 2
# Auth: u.niazi@imperial.ac.uk
# Date: 03/02/2016

source('Header.R')

#### data loading
## the dataset created earlier in previous script
load(file='Objects/lumi.n.experiment2_HJ_Habibi_040913_Sep13.rds')

# dfSymptoms = read.csv('Data_external/Original/RSV_Challenge_for_Umar/RSV_challenge_annotation_MH.csv', header=T)
# rownames(dfSymptoms) = dfSymptoms$samplesNum
# 
# dfSamples = pData(oExp.lumi)
# dfSymptoms = dfSymptoms[rownames(dfSamples),]
# i = grep('NA', rownames(dfSymptoms))
# 
# dfSamples = dfSamples[-i,]
# dfSymptoms = dfSymptoms[-i,]
# 
# # remove missing samples
# oExp.lumi = oExp.lumi[,-i]
# dfSamples$Symptoms = dfSymptoms$group
# 
# pData(oExp.lumi) = dfSamples
# rm(dfSamples)
# rm(dfSymptoms)
### data formatting done

## samples
xtabs(~ oExp.lumi$group + oExp.lumi$day)

## choose appropriate factor
## the analysis is performed along three different groups based on symptoms
## subset the data accordingly

i = grepl('^Uninfected$', x = oExp.lumi$group)
g = grepl('^0$', oExp.lumi$day)
i = g | i
## OR
i = grepl('^Asymptomatic, infected$', x = oExp.lumi$Symptoms)
g = grepl('^0$', oExp.lumi$Day)
i = g | i
## OR
i = grepl('^Clinical cold$', x = oExp.lumi$group)
g = grepl('^0$', oExp.lumi$day)
i = g | i

# select data matrix and remove genes with low detection call
ivDetection = detectionCall(oExp.lumi)
mDat = exprs(oExp.lumi)
mDat = mDat[ivDetection > 0, i]
# sample annotation
dfSamples = pData(oExp.lumi)[i,]

# sanity check
identical(dfSamples$sampleID, colnames(mDat))
table(dfSamples[ , 'group'], dfSamples[  , 'day'])

## create grouping factors by data
fSamples = rep(NA, times=nrow(dfSamples))

i = grepl('^0$', dfSamples$day)
# sanity check
table(dfSamples[i , 'group'], dfSamples[i  , 'day'])

fSamples[i] = 'D0'

# pool adjacent days together
i1 = grepl('^1$|^3$', dfSamples$day)
table(dfSamples[(i1), 'group'], dfSamples[(i1), 'day'])
i2 = grepl('^14$', dfSamples$day)
table(dfSamples[(i2), 'group'], dfSamples[(i2), 'day'])
i3 = grepl('^7$|^10$', dfSamples$day)
table(dfSamples[(i3), 'group'], dfSamples[(i3), 'day'])
i4 = grepl('^28$', dfSamples$day)
table(dfSamples[(i4), 'group'], dfSamples[(i4), 'day'])

fSamples[i1] = 'D1.3'
fSamples[i2] = 'D14'
fSamples[i3] = 'D7.10'
fSamples[i4] = 'D28'

table(fSamples)
table(is.na(fSamples))
# sanity check
data.frame(fSamples, dfSamples$day)
fSamples = factor(fSamples, levels=c('D0', 'D1.3', 'D7.10', 'D14', 'D28'))
levels(fSamples)
table(fSamples)
##### main analysis section
# map the nuID to human symbols
cvSym = getSYMBOL(rownames(mDat), 'lumiHumanAll.db')
# number of genes not annotated
table(is.na(cvSym))
# remove unannotated genes
mDat = mDat[!is.na(cvSym),]

design = model.matrix(~fSamples)
colnames(design) = levels(fSamples)
head(design)

fit = lmFit(mDat, design)
fit = eBayes(fit)

# get enterez gene id
# cvEnterez = getEG(rownames(mDat), 'lumiHumanAll.db')
# get annotation
df = select(lumiHumanAll.db, keys = rownames(mDat), columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# sanity check
nrow(mDat) == nrow(df)
# add annotation to limma object
fit$genes = df
topTable(fit, adjust='BH')

# look at top tables for each comparison
for (i in 2:length(levels(fSamples))){
  print(paste(levels(fSamples)[i], i))
  print(topTable(fit, coef=i, adjust='BH'))
}

# get the list of genes for each comparison i.e. each coefficient compared to base line
lSigGenes.adj = vector('list', length = length(levels(fSamples))-1)
names(lSigGenes.adj) = levels(fSamples)[2:length(levels(fSamples))]

for (i in 2:length(levels(fSamples))){
  p.adj = p.adjust(fit$p.value[,i], method = 'BH')
  lSigGenes.adj[[i-1]] = names(p.adj)[p.adj < 0.1]
}

#cvSigGenes.adj = unique(cvSigGenes.adj)
sapply(lSigGenes.adj, length)

######### Volcano plots
# plot volcano plots
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  f_plotVolcano(dfGenes, paste(names(n[i])), fc.lim = c(-1.5, 1.5))
}

###### HeatMaps
# make heatmaps of top genes
dfRes = topTable(fit, adjust='BH', number=Inf, p.value=0.1)
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

# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mCommonGenes.grp = cbind(mCommonGenes.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mCommonGenes.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)
# ## choose groups where there are at least 10 genes
# temp3 = temp2[which(temp2[,9] >= 10),]
# iSelGroups = temp3[,'cp']

# select groups of choice
#iSelGroups = c(1, 4, 21, 39, 45, 65, 66, 67, 70, 73, 74, 75, 77, 78, 79)
rn = rownames(mCommonGenes.grp)
#m1 = mDat[rn,]
df.rn = select(lumiHumanAll.db, keys = rn, columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# write csv to look at gene list
write.csv(df.rn[,-1], file=paste('Temp/', 'commongenes', '.csv', sep=''))

## choose appropriate combination
# ## common between 
# i = which(mCommonGenes[,'NoCold0.1'] & mCommonGenes[,'Cold0.1'])
# i = which(mCommonGenes[,'NoCold0.1'] | mCommonGenes[,'Cold0.1'])
# i = which(mCommonGenes[,'NoCold0.1'] & !mCommonGenes[,'Cold0.1'])

## all overexpressed genes if interested in
i = 1:nrow(mCommonGenes)

m1 = as.matrix(mCommonGenes[i,])
m1 = mDat[rownames(m1),]

# # if using part of data
# m1 = m1[,which(fSamples %in% c('UINC', 'NoCold0.1', 'Cold0.1'))]
# fGroups = as.character(fSamples[which(fSamples %in% c('UINC', 'NoCold0.1', 'Cold0.1'))])
# colnames(m1) = fGroups
# fGroups = factor(fGroups, levels = c('UINC', 'NoCold0.1', 'Cold0.1'))
# m1 = m1[,order(fGroups)]
# fGroups = fGroups[order(fGroups)]

# if using all data
fGroups = fSamples
colnames(m1) = fGroups
m1 = m1[,order(fGroups)]
fGroups = fGroups[order(fGroups)]

# ignore step if stabalization not required
m1 = t(apply(m1, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(m1) = fGroups

# scale across rows
rownames(m1) = dfRes[rownames(m1), 'SYMBOL']
mCounts = t(m1)
mCounts = scale(mCounts)
mCounts = t(mCounts)
# threshhold the values
mCounts[mCounts < -3] = -3
mCounts[mCounts > 3] = 3

# draw the heatmap  color='-RdBu:50'
aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)


# write results csv files
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  dfGenes.2 = dfGenes[lSigGenes.adj[[names(n[i])]],]
  rownames(dfGenes.2) = NULL
  f = paste('Results/', 'Significant_genes_at_10pcFDR_ClCold', names(n[i]), '.csv', sep='')
  dfGenes.2 = dfGenes.2[,c(2, 3, 4, 5, 6, 8, 9)]
  write.csv(dfGenes.2, file=f)
}


########### pathway analysis using CGraph library
# uniprot annotation for data
#dfGenes = topTable(fit, number = Inf)
#cvGenes = rownames(dfGenes[(dfGenes$adj.P.Val < quantile(dfGenes$adj.P.Val, probs = 0.25)),])
cvGenes = rownames(mCommonGenes)
dfGenes = select(lumiHumanAll.db, keys = cvGenes, columns = c('ENTREZID', 'SYMBOL', 'GENENAME', 'UNIPROT'), keytype = 'PROBEID')
dfGenes = na.omit(dfGenes)
#dfGenes = dfGenes[!duplicated(dfGenes$ENTREZID),]

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
# separate the factor and the count matrix
#rownames(dfGenes) = dfGenes$PROBEID
mCounts = mDat[unique(dfGenes$PROBEID),]
# map the probe names to enterez ids
i = match(rownames(mCounts), dfGenes$PROBEID)
rownames(mCounts) = dfGenes$ENTREZID[i]
fGroups = as.character(fSamples)
# select only the groups with significant genes
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
i = which(fSamples %in% levels(fSamples)[c(1, n)])
fGroups = factor(fGroups[i], levels = levels(fSamples)[c(1, n)])

# # merge the groups together
# fGroups.2 = rep(NA, length(fGroups))
# i = grep('UINC', fGroups)
# fGroups[i]
# fGroups.2[i] = 'UINC'
# i = grep('^NoCold', fGroups)
# fGroups[i]
# fGroups.2[i] = 'NoCold'
# i = grep('^Cold', fGroups)
# fGroups[i]
# fGroups.2[i] = 'Cold'
# fGroups = factor(fGroups.2, levels = c('UINC', 'NoCold', 'Cold'))

# subset the count matrix
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
i = which(fSamples %in% levels(fSamples)[c(1, n)])
mCounts = mCounts[,i]
colnames(mCounts) = fGroups
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

# stabalize the data and check correlation again
mCounts.bk = mCounts
# stabalize the data
mCounts.st = apply(mCounts, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(mCounts.st) = fGroups

# create a correlation matrix
mCor = cor(mCounts.st)
# check distribution 
hist(mCor, prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.6, bSuppressPlots = F)

## general graph structure
set.seed(1)
plot.final.graph(oGr)

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

## centrality diagnostics
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
par(p.old)

# these diagnostics plots should be looked at in combination with the centrality graphs
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)

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

dir.create('Results', showWarnings = F)
write.csv(dfTopGenes.cent, file='Results/Top_Centrality_Genes.csv')

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


## we can look at the problem from the other direction and look at clusters instead of genes
# sample plots
# mean expression of groups in every cluster
par(p.old)
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster', cex.axis=0.7)
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)
# principal component plots
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)
# plot variance of cluster
m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters
m = getClusterMarginal(oGr, t(mCounts))

csClust = rownames(m)
length(csClust)
i = 1
plot.cluster.variance(oGr, m[csClust[i:(i+1)],], fGroups, log = FALSE); i = i+2

i = 1
temp = t(as.matrix(m[csClust[i],]))
rownames(temp) = csClust[i]
plot.cluster.variance(oGr, temp, fGroups, log=FALSE); i = i+1

pdf('Temp/var.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)
dev.off(dev.cur())

i = which(dfReactome.sub$V2 %in% csClust)
dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
rownames(dfCluster.name) = NULL
dfCluster.name

# Various plots for one cluster of choice
csClust = '372790'

lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=60)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.7, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}

# heatmap of the genes
ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
mC = t(mCounts)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
mC = t(scale(t(mC)))
# threshhold the values
mC[mC < -3] = -3
mC[mC > +3] = +3
# draw the heatmap
hc = hclust(dist(mC))
aheatmap(mC, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
         annColors=NA, Colv=NA)


# if we want to plot variance of one gene at a time
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
mC = t(mCounts)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
rn = rownames(mC)
length(rn)
i = 1
pdf('Temp/var2.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, (mC), fGroups, iDrawCount = length(rn))
dev.off(dev.cur())
# temp = t(as.matrix(mC[rn[i],]))
# rownames(temp) = rn[i]
# plot.cluster.variance(oGr, temp, fGroups, log=FALSE); i = i+1


### plot all clusters in graph
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
sort(table(dfCluster$cluster))
#csClust = rownames(m$clusters)
csClust = as.character(unique(dfCluster$cluster))


# plot the graphs at each contrast
lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=30)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.14, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}
