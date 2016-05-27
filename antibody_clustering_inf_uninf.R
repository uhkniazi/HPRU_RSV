# Name: antibody_clustering_inf_uninf.R
# Auth: u.niazi@imperial.ac.uk
# Date: 03/05/2016
# Desc: import the antibody titers and cluster on fold changes


## import dataset
dfData = read.csv('Data_external/Serology data/antibody.csv', header=T)
rownames(dfData) = dfData$PID

# choose only nasal day 28 fold changes
mData = as.matrix(dfData[,4:5])

## distribution of the data
apply(mData, 2, hist)
apply(mData, 2, function(x) plot(density(x)))
summary(mData)

## cluster the data based on individual patients
hc = hclust(dist(scale(mData)))
plot(hc, main='Patient Clustering', sub='')

## try clustering with PCA
mData.pc = scale(mData)

pr.out = prcomp(mData.pc, scale=F)
biplot(pr.out)
# col.p = rainbow(nlevels(fOutcome))
# col = col.p[as.numeric(fOutcome)]
# 
plot(pr.out$x[,1:2], pch=19, xlab='Z1', ylab='Z2',
      main='PCA comp 1 and 2')

plot(mData, pch=20)
## find 3 clusters in the data
set.seed(123)
km.out = kmeans(mData, centers = 3, nstart = 20)
table(km.out$cluster)
fGroups = paste0('AB', km.out$cluster)
fGroups = factor(fGroups)
dfData$fGroups = fGroups

## boxplots for these groups
p.old = par(mfrow=c(1,2))
sapply(1:2, function(x){
  boxplot(mData[,x] ~ fGroups, main=colnames(mData)[x])
 # print(t.test(mData[,x] ~ fGroups))
})
par(p.old)

## plot all the data together, z-scaled
mData.sc = scale(mData)
# md = apply(mData, 2, mad)
# m = apply(mData, 2, median)
# 
# mData.sc = sweep(mData, 2, m, FUN = '-')
# mData.sc = sweep(mData.sc, 2, md, '/')

names(fGroups) = rownames(mData)
fGroups.or = fGroups[order(fGroups)]
mData.sc = mData.sc[order(fGroups),]

matplot(mData.sc, type='l', lty=1, xaxt='n', xlab='', ylab='Scaled Data', 
        main='Antibody titers fold change at Day 28')

axis(1, 1:nrow(mData.sc), labels = paste0(rownames(mData.sc), fGroups.or), las=2, cex.axis=0.8)

legend('topright', legend = c('Nasal.RSV', 'Nasal.F'), lty=1, col=1:2)
## test with LDA
library(MASS)
dfFit = dfData[,4:6]
lda.fit = lda(fGroups ~ ., data=dfFit)
plot(lda.fit)

write.csv(dfData, file='Results/antibody_grouping.csv')
##################### 

##################################### second subgrouping based on outcome and symptoms
dfOutcome = read.csv('Data_external/Serology data/outcome.csv', header=T, stringsAsFactors=F)
# remove some duplicates
f = duplicated(as.character(dfOutcome$PID))
dfOutcome = dfOutcome[!f,]
rownames(dfOutcome) = dfOutcome$PID

i = match(dfOutcome$PID, dfData$PID)
dfData = dfData[i,]
dfData$Outcome = dfOutcome$Outcome
dfData$Symptoms = dfOutcome$Symptoms
# 
# dfOutcome.sub = dfOutcome[rownames(dfData),]
# dfData$Outcome = dfOutcome.sub$Outcome
# dfData$Symptoms = dfOutcome.sub$Symptoms
## 

### examine expressions
## create a grouping factor on AB1, AB2, AB3 and subgrouping by outcome and symptoms
fGroups = factor(paste(dfData$fGroups, dfData$Outcome, dfData$Symptoms, sep='.'))
as.data.frame(table(fGroups))
mData = as.matrix(dfData[,4:5])

sapply(1:2, function(x){
  stripchart(mData[,x] ~ fGroups, main=colnames(mData)[x], las=2, cex.axis=0.6, vertical=T, method='jitter',
             ylab='')
  print(pairwise.t.test(mData[,x], g = fGroups, p.adjust.method = 'bonf'))
})

mData.sc = scale(mData)

names(fGroups) = rownames(mData)
fGroups.or = fGroups[order(fGroups)]
mData.sc = mData.sc[order(fGroups),]
p.old = par(mar=c(6.5, 3, 1, 1))
matplot(mData.sc[,1:2], type='l', lty=1, xaxt='n', xlab='', ylab='Scaled Data', 
        main='Antibody titers fold change at Day 28')
axis(1, 1:nrow(mData.sc), labels = paste0(rownames(mData.sc), fGroups.or), las=2, cex.axis=0.6)
legend('topright', legend = c('Nasal.RSV', 'Nasal.F'), lty=1, col=1:2)
par(p.old)

write.csv(dfData, file='Results/antibody_grouping_outcome.csv')

# 
# 
# 
# dfData.mod = data.frame(Year=factor(dfData$Year), Nasal.RSV=dfData$Nasal.RSV.IgA.fold.change.d28, 
#                         Nasal.F=dfData$Nasal.F.IgA.fold.change.d28, fGroups=fGroups)







