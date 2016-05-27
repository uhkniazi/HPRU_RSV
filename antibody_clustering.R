# Name: antibody_clustering.R
# Auth: u.niazi@imperial.ac.uk
# Date: 03/05/2016
# Desc: import the antibody titers and cluster on fold changes


## import dataset
dfData = read.csv('Data_external/Serology data/antibody.csv', header=T)
rownames(dfData) = dfData$PID

mData = as.matrix(dfData[,3:5])

## distribution of the data
apply(mData, 2, hist)
apply(mData, 2, function(x) plot(density(x)))
summary(mData)

## cluster the data based on individual patients
hc = hclust(dist(scale(mData)))
plot(hc, main='Patient Clustering', sub='')

## cluster the other direction, i.e. measurement types
hc = hclust(dist(t(scale(mData))))
plot(hc)

## try clustering with PCA
mData.pc = scale(mData)

pr.out = prcomp(mData.pc, scale=F)
biplot(pr.out)

## find 2 clusters in the data
set.seed(123)
km.out = kmeans(mData, centers = 2, nstart = 20)
table(km.out$cluster)
fGroups = paste0('AB', km.out$cluster)
fGroups = factor(fGroups)
dfData$fGroups = fGroups

## boxplots for these groups
sapply(1:3, function(x){
  boxplot(mData[,x] ~ fGroups, main=colnames(mData)[x])
  print(t.test(mData[,x] ~ fGroups))
})

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

legend('topright', legend = c('Serum', 'Nasal.RSV', 'Nasal.F'), lty=1, col=1:3)
## test with LDA
library(MASS)
dfFit = dfData[,3:6]
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

dfOutcome.sub = dfOutcome[rownames(dfData),]
dfData$Outcome = dfOutcome.sub$Outcome
dfData$Symptoms = dfOutcome.sub$Symptoms
## 

### examine expressions
fGroups = factor(paste0(dfData$fGroups, dfData$Symptoms))
mData = as.matrix(dfData[,3:5])

sapply(1:3, function(x){
  boxplot(mData[,x] ~ fGroups, main=colnames(mData)[x])
  print(pairwise.t.test(mData[,x], g = fGroups, p.adjust.method = 'bonf'))
})

mData.sc = scale(mData)

names(fGroups) = rownames(mData)
fGroups.or = fGroups[order(fGroups)]
mData.sc = mData.sc[order(fGroups),]
matplot(mData.sc[,2:3], type='l', lty=1, xaxt='n', xlab='', ylab='Scaled Data', 
        main='Antibody titers fold change at Day 28')
axis(1, 1:nrow(mData.sc), labels = paste0(rownames(mData.sc), fGroups.or), las=2, cex.axis=0.8)
legend('topright', legend = c('Nasal.RSV', 'Nasal.F'), lty=1, col=1:2)







dfData.mod = data.frame(Year=factor(dfData$Year), Nasal.RSV=dfData$Nasal.RSV.IgA.fold.change.d28, 
                        Nasal.F=dfData$Nasal.F.IgA.fold.change.d28, Symptoms=factor(dfData$Symptoms), 
                        fGroups=dfData$fGroups)







