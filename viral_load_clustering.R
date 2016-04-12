# Name: viral_load_clustering.R
# Auth: u.niazi@imperial.ac.uk
# Date: 12/04/2016
# Desc: import of viral load dataset, and grouping patients based on expression levels


## import dataset
dfData = read.csv('Data_external/RSV 2011-13 infected VL.csv', header=T)
rownames(dfData) = dfData$Days
dfData = dfData[,-1]

# paremeters from the data
mData = as.matrix(dfData)
days = rowMeans(mData, na.rm = T)
fDays = cut(as.integer(names(days)), breaks = 5, include.lowest = T, labels = 1:5) 
table(fDays)
meanDays = tapply(days, fDays, mean, na.rm=T)
nDays = length(days)
Sdays = var(days, na.rm=T) * (nDays-1)

# assign factor to table
dfData$fDays = fDays

## calculate new data for each patient
get.post.pred = function(dat, param){
  S = param$var
  n = param$length
  # sample variance
  sig2 = S/rchisq(1, n-1)
  dat = c(dat, param$binMean)
  # sample mean
  mu = rnorm(1, mean(dat, na.rm=T), sqrt(sig2)/sqrt(n))
  #return(rnorm(1, mu, sqrt(sig2)))
  return(mu)
}

lparam = list(var=Sdays, length=nDays)

mData.new = apply(mData, 2, function(x){
  y = NULL
  for (i in 1:nlevels(dfData$fDays)){
    param = lparam
    param$binMean = meanDays[i]
    y = c(y, get.post.pred(x[dfData$fDays == levels(dfData$fDays)[i]], param))
  }
  return(y)
})

hc = hclust(dist(scale(mData.new)))
plot(hc)

hc = hclust(dist(t(scale(mData.new))))
plot(hc)

mData.pc = scale(mData.new)

pr.out = prcomp(mData.pc, scale=F)
biplot(pr.out)

plot(hclust(dist(mData.pc)))
matplot(mData.pc, type='l', xlab='bins', ylab='expressions')

mData.pc = t(mData.pc)
pr.out = prcomp(mData.pc, scale=F)
biplot(pr.out)

plot(hclust(dist(mData.pc)))



tapply(dfData[,1], dfData$fDays, get.post.pred, lparam)




