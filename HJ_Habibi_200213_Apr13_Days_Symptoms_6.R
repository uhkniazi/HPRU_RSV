# File: HJ_Habibi_200213_Apr13_Days_Symptoms_6.R
# Desc: Data set from RSV challenge - using days and symptoms as grouping factors 
# Auth: u.niazi@imperial.ac.uk
# Date: 2/11/15

source('Header.R')

#### data loading
## the dataset created earlier in previous script
load(file='Objects/lumi.n.experiment1_HJ_Habibi_200213_Apr13.rds')

dfSymptoms = read.csv('Data_external/Original/RSV_Challenge_for_Umar/RSV_challenge_annotation_MH.csv', header=T)
rownames(dfSymptoms) = dfSymptoms$samplesNum

dfSamples = pData(oExp.lumi)
dfSymptoms = dfSymptoms[rownames(dfSamples),]
i = grep('NA', rownames(dfSymptoms))

dfSamples = dfSamples[-i,]
dfSymptoms = dfSymptoms[-i,]

# remove missing samples
oExp.lumi = oExp.lumi[,-i]
dfSamples$Symptoms = dfSymptoms$group
# remove day 0
i = which(dfSamples$Day == '0')
oExp.lumi = oExp.lumi[,-i]
dfSamples = dfSamples[-i,]

pData(oExp.lumi) = dfSamples

xtabs(~ dfSamples$Symptoms + dfSamples$Day)

## choose appropriate factor
fSamples = rep(NA, length=nrow(pData(oExp.lumi)))

i = grepl('^Uninfected$', x = oExp.lumi$Symptoms)
# group days 0 uninfected group
# i2 = grepl('^0$', x = oExp.lumi$Day)
# fSamples[i | i2 ] = 'UINC'

fSamples[i] = 'UINC'

# sanity check
table(dfSamples[i , 'Symptoms'], dfSamples[i  , 'Day'])

# create second level by grouping 
# infected together merged and split into 3 days
i = grepl('^Clinical cold$', x = oExp.lumi$Symptoms)
# pool adjacent days together
i1 = grepl('^1$', x = oExp.lumi$Day)
i2 = grepl('^3$|^5$', x = oExp.lumi$Day)
i3 = grepl('^7$|^10$', x = oExp.lumi$Day)
i4 = grepl('^28$', x = oExp.lumi$Day)

# sanity check
table(dfSamples[i & (i1 | i2 | i3 | i4), 'Symptoms'], dfSamples[i & (i1 | i2 | i3 | i4), 'Day'])
fSamples[i & i1] = paste0('Cold', '1')
fSamples[i & i2] = paste0('Cold', '3.5')
fSamples[i & i3] = paste0('Cold', '7.10')
fSamples[i & i4] = paste0('Cold', dfSamples[i & i4, 'Day'])

# infected no cold together merged and split into 4 days
i = grepl('^Asymptomatic, infected$', x = oExp.lumi$Symptoms)
# pool adjacent days together
i1 = grepl('^1$', x = oExp.lumi$Day)
i2 = grepl('^3$|^5$', x = oExp.lumi$Day)
i3 = grepl('^7$|^10$', x = oExp.lumi$Day)
i4 = grepl('^28$', x = oExp.lumi$Day)

# sanity check
table(dfSamples[i & (i1 | i2 | i3 | i4), 'Symptoms'], dfSamples[i & (i1| i2 | i3 | i4), 'Day'])
fSamples[i & i1] = paste0('NoCold', '1')
fSamples[i & i2] = paste0('NoCold', '3.5')
fSamples[i & i3] = paste0('NoCold', '7.10')
fSamples[i & i4] = paste0('NoCold', dfSamples[i & i4, 'Day'])

table(fSamples)

fSamples = factor(fSamples, levels = c('UINC', 'NoCold1', 'NoCold3.5', 'NoCold7.10', 'NoCold28', 
                                      'Cold1', 'Cold3.5', 'Cold7.10', 'Cold28'))
as.data.frame(table(fSamples))
############# Select genes based on differential expression analysis.
mDat = exprs(oExp.lumi)

# select genes with low detection call
ivDetection = detectionCall(oExp.lumi)
mDat = mDat[ivDetection > 0,]
# map the nuID to human symbols
cvSym = getSYMBOL(rownames(mDat), 'lumiHumanAll.db')
# number of genes not annotated
table(is.na(cvSym))
# remove unannotated genes
mDat = mDat[!is.na(cvSym),]

design = model.matrix(~fSamples)
colnames(design) = levels(fSamples)

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
  #dfGenes.2 = dfGenes[lSigGenes.adj[[names(n[i])]],]
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < 0.1)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=paste(names(n[i])), xlim=c(-3, 3))
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
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
# create groups in the data based on 2^8-1 combinations
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

## choose groups where there are at least 10 genes
temp3 = temp2[which(temp2[,9] >= 10),]
iSelGroups = temp3[,'cp']

# select groups of choice
#iSelGroups = c(1, 4, 21, 39, 45, 65, 66, 67, 70, 73, 74, 75, 77, 78, 79)
rn = rownames(mCommonGenes.grp[cp %in% iSelGroups[c(1, 3, 4)],])
m1 = mDat[rn,]
df.rn = select(lumiHumanAll.db, keys = rn, columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# write csv to look at gene list
write.csv(df.rn[,-1], file=paste('Temp/', 'NoCold3.5and7.10', '.csv', sep=''))

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
  f = paste('Results/', 'Significant_genes_at_10pcFDR_NoD0_UINC_vs_', names(n[i]), '.csv', sep='')
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
dfGenes = dfGenes[!duplicated(dfGenes$ENTREZID),]

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
rownames(dfGenes) = dfGenes$PROBEID
mCounts = mDat[rownames(dfGenes),]
rownames(mCounts) = dfGenes$ENTREZID
fGroups = as.character(fSamples)
# select only the groups with significant genes
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
i = which(fSamples %in% levels(fSamples)[c(1, n)])
fGroups = factor(fGroups[i], levels = levels(fSamples)[c(1, n)])
# subset the count matrix
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
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# stabalize the data and check correlation again
mCounts.bk = mCounts
# stabalize the data
mCounts = apply(mCounts, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(mCounts) = fGroups

# create a correlation matrix
mCor = cor(mCounts)
# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.5, bSuppressPlots = F)

## general graph structure
set.seed(1)
plot.final.graph(oGr)

## community structure
# plot the main communities and centrality graphs
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 20)
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
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 40)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)

# these diagnostics plots should be looked at in combination with the centrality graphs
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)

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

write.csv(dfTopGenes.cent, file='Results/Top_Centrality_Genes_HJ_Habibi_200213_Apr13_NC_vs_Cold.csv')

# plot a heatmap of these top genes
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


# plot a graph of these top genes
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = as.character(dfTopGenes.cent$VertexID))
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize=30)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.2, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', main='Cold7.10 vs NC')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

# switch the factor levels
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = as.character(dfTopGenes.cent$VertexID))
fG = factor(fGroups, levels = c('NC', 'Cold7.10', 'Cold3.5'))
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fG, bColor = T, iSize=30)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.2, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', main='Cold3.5 vs NC')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))




## location of the largest clique
# plot largest connected graph - clique
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 100)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

### largest clique 
# plot the graph with clique highlighted
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 40)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, layout=layout_with_fr)
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

# repeat after changing levels
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fG, bColor = T, iSize = 40)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, layout=layout_with_fr)
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
par(p.old)


## we can look at the problem from the other direction and look at clusters instead of genes
# sample plots
# mean expression of groups in every cluster
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'topright', main='Total Change in Each Cluster')
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)

# plot variance of cluster
m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters
plot.cluster.variance(oGr, m[c('1280218', '1280215'),], fGroups)

plot.cluster.variance(oGr, m[csClust[i:(i+1)],], fGroups); i = i+2


#### plot a graph of top clusters clusters 
m = getSignificantClusters(oGr, t(mCounts), fGroups)
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
sort(table(dfCluster$cluster))
csClust = names(sort(table(dfCluster$cluster), decreasing = T))
csClust = csClust[csClust %in% rownames(m)]
#csClust = rownames(m$clusters)

# plot these genes 
# plot for each factor level
lev = levels(fGroups)[-1]
m = mCounts
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=50)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.14, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
       main=paste(lev[i], 'vs UINC'))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}

dfCluster = dfCluster[dfCluster$cluster %in% csClust,]
df = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
dfCluster = cbind(dfCluster[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
write.csv(dfCluster, file='Results/Top_Clusters_HJ_Habibi_200213_Apr13_UINC_vs_SymptomsAndDays.csv')

### Graphs of clusters
# get clusters of choice to make subplots
sort(table(dfCluster$cluster))

# plot one cluster of choice
# plot one cluster of choice
csClust = c('1280218')
set.seed(1)
ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust)
#ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), fGroups, bColor = T, iSize = 50)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
# V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
# plot(ig.sub, vertex.label.cex=0.7, layout=layout_with_fr)
# # switch grouping levels
# ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), fG, bColor = T, iSize = 50)
# set.seed(1)
# plot(ig.sub, vertex.label.cex=0.7, layout=layout_with_fr)

mC = t(mCounts)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
mC = t(mC)
mC = apply(mC, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(mC) = rownames(mCounts)
mC = t(mC)
mC = t(scale(t(mC)))
# threshhold the values
mC[mC < -3] = -3
mC[mC > +3] = +3
# draw the heatmap
hc = hclust(dist(mC))
aheatmap(mC, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
         annColors=NA, Colv=NA)



