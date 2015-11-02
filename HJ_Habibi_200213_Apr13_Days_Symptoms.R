# File: HJ_Habibi_200213_Apr13_Days_Symptoms.R
# Desc: Data set from RSV challenge - using days and symptoms as grouping factors 
# Auth: u.niazi@imperial.ac.uk
# Date: 2/11/15

source('Header.R')

#### data loading
## the dataset created earlier in previous script
load(file='Objects/lumi.n.experiment1_HJ_Habibi_200213_Apr13.rds')

dfSymptoms = read.csv('Data_external/Original/RSV_Challenge_for_Umar/RSV_challange_annotation_v3.csv', header=T)
rownames(dfSymptoms) = dfSymptoms$Samples

dfSamples = pData(oExp.lumi)
dfSymptoms = dfSymptoms[rownames(dfSamples),]
i = grep('NA', rownames(dfSymptoms))

dfSamples = dfSamples[-i,]
dfSymptoms = dfSymptoms[-i,]

# remove missing samples
oExp.lumi = oExp.lumi[,-i]
dfSamples$Symptoms = dfSymptoms$Symptoms
pData(oExp.lumi) = dfSamples

## choose appropriate factor
fSamples = rep(NA, length=nrow(pData(oExp.lumi)))

i = grepl('Uninfected no cold', x = oExp.lumi$Symptoms)
# group days 0 and 1 into uninfected group
i2 = grepl('^0$|^1$', x = oExp.lumi$Day)
fSamples[i | i2 ] = 'UINC'

# sanity check
table(dfSamples[i | i2, 'Symptoms'], dfSamples[i | i2, 'Day'])

# create second level by grouping 
# infected together merged and split into 3 days
i = !grepl('Uninfected no cold', x = oExp.lumi$Symptoms)
# these days have been grouped into uninfected group
i2 = !grepl('^0$|^1$', x = oExp.lumi$Day)

# sanity check
table(dfSamples[i & i2, 'Symptoms'], dfSamples[i & i2, 'Day'])
fSamples[i & i2] = paste0('INF', dfSamples[i & i2, 'Day'])
table(fSamples)

fSamples = factor(fSamples, levels = c('UINC', 'INF3', 'INF5', 'INF7', 'INF10', 'INF28'))
table(fSamples)
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
n = (which(sapply(lSigGenes.adj, length) > 10)) + 1

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
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=paste(names(n[i])), xlim=c(-2, 2))
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
n = (which(sapply(lSigGenes.adj, length) > 10)) + 1
cvCommonGenes = NULL
for (i in seq_along(n)) {
  cvCommonGenes = append(cvCommonGenes, lSigGenes.adj[[names(n[i])]])
}
cvCommonGenes = unique(cvCommonGenes)
mCommonGenes = sapply(seq_along(n), function(x) cvCommonGenes %in% lSigGenes.adj[[names(n[x])]])
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = names(n)

## choose appropriate combination
## genes common amongs all stages
#i = which(rowSums(mCommonGenes) == 2)
## all overexpressed genes
i = 1:nrow(mCommonGenes)

m1 = as.matrix(mCommonGenes[i,])
m1 = mDat[rownames(m1),]
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
n = (which(sapply(lSigGenes.adj, length) > 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  dfGenes.2 = dfGenes[lSigGenes.adj[[names(n[i])]],]
  rownames(dfGenes.2) = NULL
  f = paste('Results/', 'Significant_genes_at_10pcFDR_UINC_vs_', names(n[i]), '.csv', sep='')
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
fGroups = fSamples
colnames(mCounts) = fGroups
mCounts = mCounts[,order(fGroups)]
fGroups = fGroups[order(fGroups)]
mCounts = t(mCounts)

# select genes that have a reactome term attached
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
plot.centrality.diagnostics(oGr)

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

write.csv(dfTopGenes.cent, file='Results/Top_Centrality_Genes_HJ_Habibi_200213_Apr13_UINC_vs_INF.csv')

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
ig = induced_subgraph(getFinalGraph(oGr), vids = dfTopGenes.cent$VertexID)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize=100)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.3, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', main='IC vs UINC')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

# switch the factor levels
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = dfTopGenes.cent$VertexID)
fG = factor(fGroups, levels = c('UINC', 'IC', 'INC'))
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fG, bColor = T, iSize=100)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.3, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', main='INC vs UINC')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))




## location of the largest clique
# plot largest connected graph - clique
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 100)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, layout=layout_with_fr)
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

# repeat after changing levels
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fG, bColor = T)
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
par(mar=c(7, 3, 2, 2)+0.1)
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomright', main='Total Change in Each Cluster')
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=2, bStabalize = F, cex.axis=0.7)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = F)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# expression in all clusters
plot.heatmap.all(oGr, t(mCounts))
# marginal expression level in each cluster
plot.heatmap.marginal(oGr, t(mCounts))
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)

dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
dfCluster = dfCluster[order(dfCluster$cluster),]
write.csv(dfCluster, file='Results/Top_Clusters_HJ_Habibi_200213_Apr13_Symptoms.csv')

### Graphs of clusters
# get clusters of choice to make subplots
sort(table(dfCluster$cluster))

# plot one cluster of choice
csClust = '372790'
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust)
ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), fGroups, bColor = T)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
plot(ig.sub, vertex.label.cex=0.7, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', main='IC vs UINC')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

# switch grouping levels
ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), fG, bColor = T)
set.seed(1)
plot(ig.sub, vertex.label.cex=0.7, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', main='INC vs UINC')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

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


