x.uninf = uninfected[[3]][,'4600']
f.uninf = uninfected[[2]]
x.cold = cold[[3]][,'4600']
f.cold = cold[[2]]
boxplot(x.uninf ~ f.uninf)
boxplot(x.cold ~ f.cold)

source('../CCrossValidation/CCrossValidation.R')
fGroups.bk = fGroups
fGroups = as.character(fGroups)
i = which(fGroups == 'UI7.10')
fGroups[i] = 'D0'
i = which(fGroups != 'D0')
fGroups[i] = 'Cold'

g = as.character(dfTopGenes.cent$VertexID)
mDat = mCounts[,g]

test = sample(1:length(fGroups), size = length(fGroups)*0.2, replace = F)
dfData = as.data.frame(mDat)
fGroups = factor(fGroups, levels=c('D0', 'Cold'))

oVar.r = CVariableSelection.RandomForest(dfData[-test,], fGroups[-test])
plot.var.selection(oVar.r)

dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)

cvTopGenes = rownames(dfRF)[1:30]
cvTopGenes = gsub('X(\\d+)', '\\1', cvTopGenes)

dfData.2 = dfData[,colnames(dfData) %in% cvTopGenes]
oVar.sub = CVariableSelection.ReduceModel(dfData.2, fGroups, boot.num = 10)

g1 = as.numeric(dfDat[2,2:24])
dat = data.frame(g1, fGroups)

l = f_lsimpost(g1[fGroups == 'Inf'])

rnorm(1, median(l$mu), sqrt(median(l$var)))
simpost(g1[fGroups == 'Inf'])





dfui.7 = read.csv(file.choose(), header = T, row.names=1)
dfui.10 = read.csv(file.choose(), header = T, row.names=1)
dfui.28 = read.csv(file.choose(), header = T, row.names=1)

dfnc.3.5 = read.csv(file.choose(), header = T, row.names=1)
dfnc.7.10 = read.csv(file.choose(), header = T, row.names=1)


dfReactome.terms = dfReactome.sub[!duplicated(dfReactome.sub$V2),]
dfReactome.terms[dfReactome.terms$V2 == '1280218', c('V2', 'V4')]

dfReactome.terms[dfReactome.terms$V2 %in% csClust, c('V2', 'V4')]






library(annotate)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(lumi)
library(limma)
source('Header.R')

source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')

dfSamples = read.csv(file.choose(), header=T)
rownames(dfSamples) = dfSamples$Samples

csFile = 'Data_external/Original/RSV_Challenge_for_Umar/RSV_challange_handover(Artem)/original annotation and data files/HJ_Habibi_200213_Apr13_No uRNA_Bkg sub_No Norm_Sample_Probe_Profile.txt'
x.lumi = lumiR.batch(csFile, lib.mapping = 'lumiHumanIDMapping')

# add the sample annotation
# check if sample ids same
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
fSamples = as.factor(x.lumi$Day)
fSamples = as.factor(x.lumi$Sample.group.3)
fSamples = as.factor(x.lumi$Study.Group)
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
fSamples = as.factor(lumi.n.q$Day)
fSamples = as.factor(lumi.n.q$Sample.group.3)
fSamples = as.factor(lumi.n.q$Study.Group)

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

# remove the outliers 
# remove the outlier groups from the data
# these can be seen on the pc2 and pc3 plots
m = pr.out$x[,1:3]
m = data.frame(m, fSamples)
i = which(m$PC2 > 0 & m$PC1 < -300)
i = unique(c(i, which(m$PC3 < -50)))
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
pData(lumi.n.q)[i,]

# remove the outliers
oExp.lumi = lumi.n.q[,-i]

# plot another PCA
m = exprs(oExp.lumi)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(oExp.lumi$Day)
fSamples = as.factor(oExp.lumi$Sample.group.3)
fSamples = as.factor(oExp.lumi$Study.Group)

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

# save the normalized object
save(oExp.lumi, file='Objects/lumi.n.experiment1_HJ_Habibi_200213_Apr13.rds')

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
# choose the factor
fSamples = factor(oExp.lumi$Study.Group, labels=c('t1', 't2', 't3'))
fSamples = as.factor(oExp.lumi$Day)

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
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=paste(names(n[i])), xlim=c(-1, 1))
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
## genes common amongs all three stages
i = which(rowSums(mCommonGenes) == 3)
## common between days 7 and 10
i = which(mCommonGenes[,1] & mCommonGenes[,2])
## all overexpressed genes
i = 1:nrow(mCommonGenes)

m1 = mCommonGenes[i,]
m1 = mDat[rownames(m1),]
m1 = m1[,which(fSamples %in% c('0', '7', '10', '28'))]
fGroups = as.character(fSamples[which(fSamples %in% c('0', '7', '10', '28'))])
colnames(m1) = fGroups
fGroups = factor(fGroups, levels = c('0', '7', '10', '28'))
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
  f = paste('Results/', 'Significant_genes_at_10pcFDR_day0_vs_', names(n[i]), '.csv', sep='')
  dfGenes.2 = dfGenes.2[,c(2, 3, 4, 5, 6, 8, 9)]
  write.csv(dfGenes.2, file=f)
}


########### pathway analysis using CGraph library
# uniprot annotation for data
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
mCounts = mCounts[,which(fSamples %in% c('0', '7', '10', '28'))]
fGroups = as.character(fSamples[which(fSamples %in% c('0', '7', '10', '28'))])
colnames(mCounts) = fGroups
fGroups = factor(fGroups, levels = c('0', '7', '10', '28'))
mCounts = mCounts[,order(fGroups)]
fGroups = fGroups[order(fGroups)]

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[n,]
mCounts = t(mCounts)
levels(fGroups)

# create correlation matrix
mCor = cor(mCounts)
# check distribution 
hist(sample(mCor, 10000, replace = T), prob=T, main='Correlation of genes', xlab='', family='Arial')

# create graph object manually
save(oGr, file='Objects/lumi.n.experiment1_HJ_Habibi_200213_Apr13_graph_object.rds')
# reload from here next time
# plots 
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=2, bStabalize = T, cex.axis=0.7)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)

# get clusters of choice to make subgraphs
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
dfCluster = dfCluster[order(dfCluster$cluster),]

# plot one cluster of choice on pdf
csClust = '168256'
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust)
ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), factor(as.character(fGroups)), bColor = T, iSize = 1000)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
plot(ig.sub, vertex.label.cex=0.7, layout=layout_with_fr)
par(p.old)

mC = t(mCounts)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL

# stabalize the data
mC = t(apply(mC, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(mC) = fGroups
mC = t(scale(t(mC)))
# threshhold the values
mC[mC < -3] = -3
mC[mC > +3] = +3
# draw the heatmap
hc = hclust(dist(mC))
aheatmap(mC, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
         annColors=NA, Colv=NA)

# table for printing
dfCluster = dfCluster[dfCluster$cluster == csClust,]
dfCluster = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
cvSum = f_csGetGeneSummaryFromGenbank(dfCluster$ENTREZID)
cvSum.2 = dfCluster$SYMBOL
dfCluster$Summary = cvSum[cvSum.2]
write.csv(dfCluster, file='Temp/dfCluster_immune.csv')





# p.adj = p.adjust(fit$p.value[,2], method = 'BH')
# cvSigGenes.adj = rownames(mDat)[p.adj < 0.1]
dfRes = topTable(fit, adjust='BH', number=Inf, p.value=0.1)
# do one level at a time
cvSigGenes.adj = lSigGenes.adj[['10']]
dfRes = dfRes[cvSigGenes.adj,]
dfRes = dfRes[order(dfRes$adj.P.Val, decreasing = F),]


############## go stats
library(GOstats)
cvSigGenes.adj.ent = unique(dfRes$ENTREZID)
params = new('GOHyperGParams', geneIds=cvSigGenes.adj.ent,
             annotation='lumiHumanAll.db',
             ontology='BP',
             pvalueCutoff= 0.01,
             conditional=FALSE,
             testDirection='over')

oGOStat = hyperGTest(params) 
# get pvalues
ivPGO = pvalues(oGOStat)
# fdr
ivPGO.adj = p.adjust(ivPGO, 'BH')

cvSigGO.ID = names(ivPGO.adj[ivPGO.adj < 0.01])
dfSigGO = summary(oGOStat)
dfSigGO = dfSigGO[dfSigGO$GOBPID %in% cvSigGO.ID,]
#dfSigGO = getGOTerm(cvSigGO.ID)[['BP']]






############### graph analysis
#library(reactome.db)
#library(org.Hs.eg.db)
dfMap = dfRes
source('../CGraphClust/CGraphClust.R')

### use reactome terms instead of go
dfMap = na.omit(dfRes)
url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
dir.create('Data_external', showWarnings = F)
csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
# map reactome pathways
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
x = gsub('\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfMap$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfMap$UNIPROT)
dfReactome.sub$ENTREZID = dfMap$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)
# assign reactome ids
# dfGraph = AnnotationDbi::select(reactome.db, rownames(mCounts), 'REACTOMEID', 'ENTREZID')
# dfGraph = na.omit(dfGraph)

### use go terms instead of reactome
dfGO = select(lumiHumanAll.db, unique(dfMap$ENTREZID), 'GO', keytype = 'ENTREZID' )
head(dfGO)
dfGO = dfGO[dfGO$ONTOLOGY == 'BP', ]
head(dfGO)
dfGO = dfGO[,1:2]
dfGraph = na.omit(dfGO)
head(dfGraph)

# separate the factor and the count matrix
fGroups = fSamples
mCounts = mDat[rownames(dfMap),]
rownames(mCounts) = dfMap$ENTREZID
mCounts = mCounts[!duplicated(rownames(mCounts)),]

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[n,]
mCounts = t(mCounts)
levels(fSamples)
# choose base level and another level
i = which(fSamples %in% c('0', '10'))
mCounts = mCounts[i,]
fGroups = as.factor(as.character(fSamples[i]))

# create a correlation matrix

# stabalize the data
mCounts = apply(mCounts, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(mCounts) = fGroups

# create a correlation matrix
mCor = cor(mCounts)
# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial')

# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.5)

# order the count matrix before making heatmaps or plots
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

# plot the main communities and centrality graphs
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 400)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
set.seed(1)
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 200)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')

# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)
plot.centrality.diagnostics(oGr)

# make a pdf output
pdf('Temp/Figures/Graph_structure.pdf')
par(mar=c(1,1,1,1)+0.1, family='Helvetica')
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
set.seed(1)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')
set.seed(1)
ig = plot.centrality.graph(oGr)
dev.off(dev.cur())


# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)
# top  5% of the vertices from each category and largest clique
l = lGetTopVertices(oGr, iQuantile = 0.85)
# top genes based on centrality parameters
dfClique = f_dfGetGeneAnnotation(l$clique)
cvSum = f_csGetGeneSummaryFromGenbank(dfClique$ENTREZID)
cvSum.2 = dfClique$SYMBOL
dfClique$Summary = cvSum[cvSum.2]
write.csv(dfClique, file='Temp/dfClique.csv')

dfDegree = f_dfGetGeneAnnotation(l$degree)
cvSum = f_csGetGeneSummaryFromGenbank(dfDegree$ENTREZID)
cvSum.2 = dfDegree$SYMBOL
dfDegree$Summary = cvSum[cvSum.2]
write.csv(dfDegree, file='Temp/dfDegree.csv')

dfHub = f_dfGetGeneAnnotation(l$hub)
cvSum = f_csGetGeneSummaryFromGenbank(dfHub$ENTREZID)
cvSum.2 = dfHub$SYMBOL
dfHub$Summary = cvSum[cvSum.2]
write.csv(dfHub, file='Temp/dfHub.csv')

dfBetweenness = f_dfGetGeneAnnotation(l$betweenness)
cvSum = f_csGetGeneSummaryFromGenbank(dfBetweenness$ENTREZID)
cvSum.2 = dfBetweenness$SYMBOL
dfBetweenness$Summary = cvSum[cvSum.2]
write.csv(dfBetweenness, file='Temp/dfBetweenness.csv')

dfCloseness = f_dfGetGeneAnnotation(l$closeness)
cvSum = f_csGetGeneSummaryFromGenbank(dfCloseness$ENTREZID)
cvSum.2 = dfCloseness$SYMBOL
dfCloseness$Summary = cvSum[cvSum.2]
write.csv(dfCloseness, file='Temp/dfCloseness.csv')


# plot largest connected graph - clique
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 50)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 3000)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
plot(ig, layout=layout_with_fr)
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
par(p.old)

# plot the graphs of centrality parameter genes 
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = dfDegree$ENTREZID)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.5, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey', main='Degree')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

ig = induced_subgraph(getFinalGraph(oGr), vids = dfCloseness$ENTREZID)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.5, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey', main='Closeness')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

ig = induced_subgraph(getFinalGraph(oGr), vids = dfBetweenness$ENTREZID)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.5, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey', main='Betweenness')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

ig = induced_subgraph(getFinalGraph(oGr), vids = dfHub$ENTREZID)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.5, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey', main='Hub')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
par(p.old)

## we can look at the problem from the other direction and look at clusters instead of genes
# sample plots
# mean expression of groups in every cluster
par(mar=c(7, 3, 2, 2)+0.1)
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomright', main='Total Change in Each Cluster')
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=2, bStabalize = T, cex.axis=0.7)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# expression in all clusters
plot.heatmap.all(oGr, t(mCounts))
# marginal expression level in each cluster
plot.heatmap.marginal(oGr, t(mCounts))
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)

# # plot selected clusters
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '1280215')
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '168249')
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '3247509')
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '72203')


library(GO.db)
columns(GO.db)
temp = getSignificantClusters(oGr, t(mCounts), fGroups, bStabalize = T)
rn = rownames(temp$clusters)
select(GO.db, keys = rn, columns = 'DEFINITION', keytype = 'GOID')
temp = select(GO.db, keys = rn, columns = 'DEFINITION', keytype = 'GOID')


# get clusters of choice to make subgraphs
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
dfCluster = dfCluster[order(dfCluster$cluster),]

# plot one cluster of choice on pdf
csClust = '1280215'
set.seed(1)
ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust)
ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), fGroups, bColor = T)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
plot(ig.sub, vertex.label.cex=0.7, layout=layout_with_fr)
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

# table for printing
csClust = '1280215'
dfCluster = dfCluster[dfCluster$cluster == csClust,]
dfCluster = f_dfGetGeneAnnotation(as.character(dfCluster$gene))

l = lGetTopVertices(oGr, iQuantile = 0.85)
sapply(l, length)
f = sapply(seq_along(1:5), function(x) dfCluster$ENTREZID %in% l[[x]])
colnames(f) = names(l)
dfCluster = cbind(dfCluster, f)
n = f_csGetGeneSummaryFromGenbank(dfCluster$ENTREZID)
cvSum.2 = as.character(dfCluster$SYMBOL)
dfCluster$Summary = n[cvSum.2]
write.csv(dfCluster, file='Temp/Figures/cluster_1280215.csv')


# plot the subgraphs of the significant clusters
l = getSignificantClusters(oGr, mCounts = t(mCounts), fGroups, bStabalize = T)
csClust = rownames(l$clusters)
pdf('Temp/graphs.pdf')
par(mar=c(1,1,1,1)+0.1)
sapply(seq_along(csClust), function(x){
  set.seed(1)
  ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust[x])
  ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), fGroups, bColor = T)
  n = f_dfGetGeneAnnotation(V(ig.sub)$name)
  V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
  plot(ig.sub, vertex.label.cex=0.5, layout=layout_with_fr, main=csClust[x])
  ig.sub = getLargestCliqueInCluster(oGr, csClustLabel = csClust[x])
  ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), fGroups, bColor = T)
  n = f_dfGetGeneAnnotation(V(ig.sub)$name)
  V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
  plot(ig.sub, vertex.label.cex=0.8, layout=layout_with_fr, main=csClust[x])
  plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = csClust[x])
})
sym = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
cluster = as.character(dfCluster$cluster)
dfCluster = cbind(sym, cluster)
cvSum = f_csGetGeneSummaryFromGenbank(iID = as.character(dfCluster$ENTREZID))
cvSum.2 = as.character(dfCluster$SYMBOL)
dfCluster$Summary = cvSum[cvSum.2]

write.csv(dfCluster, file='Temp/dfCluster.csv')

# saving graph object to visualize in cytoscape or other graph viewers
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 20)
set.seed(1)
par(mar=c(1,1,1,1)+0.1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color='grey')
set.seed(1)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color='grey')

n = f_dfGetGeneAnnotation(V(ig)$name)
V(ig)[n$ENTREZID]$label = n$SYMBOL
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
dir.create('Test_data', showWarnings = F)
write.csv(dfCluster, 'Test_data/clusters.csv')
write.graph(ig, file = 'Test_data/graph.graphml', format = 'graphml')



############## end graph






