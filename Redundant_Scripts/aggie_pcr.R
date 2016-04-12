# Name: aggie_pcr.R
# Date: 15/12/15
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Desc: analysis of pcr data


# data input
dfDat = read.csv('Data_external/aggie_pcr/input.csv', header=T)


fGroups = c(rep('Inf', times=14), rep('Ninf', times=9))
fGroups = factor(fGroups, levels = c('Ninf', 'Inf'))


## function to simulate posterior predictive for missing data
simpost = function(ivDat){
  #set.seed(123)
  # calculate using non-informative prior parameters  
  sigma.0 = 0
  k.0 = 0
  v.0 = k.0 - 1
  mu.0 = 0
  
  ## look at page 68 of Bayesian Data Analysis (Gelman) for formula
  sim.post = function(dat.grp){
    # calculate conjugate posterior
    n = length(dat.grp)
    k.n = k.0 + n
    v.n = v.0 + n
    y.bar = mean(dat.grp)
    s = sd(dat.grp)
    mu.n = (k.0/k.n * mu.0) + (n/k.n * y.bar)
    sigma.n = (( v.0*sigma.0 ) + ( (n-1)*(s^2) ) + ( (k.0*n/k.n)*((y.bar-mu.0)^2) )) / v.n
    #post.scale = ((prior.dof * prior.scale) + (var(dat.grp) * (length(dat.grp) - 1))) / post.dof
    ## simulate variance
    sigma = (sigma.n * v.n)/rchisq(1, v.n)
    mu = rnorm(1, mu.n, sqrt(sigma)/sqrt(k.n))
    return(list(mu=mu, var=sigma))
  }
  
  # get a sample of the posterior variance for each group
  s = sim.post(ivDat)
  return(rnorm(1, s$mu, sqrt(s$var)))
}

############ end function

# assign unique names
rn = paste(dfDat$X, 1:nrow(dfDat), sep='.')
rownames(dfDat) = rn
dfDat = dfDat[,-1]

# divide by housekeeping gene
m = as.matrix(dfDat)
h.k = m['RTC.91',]
m2 = sweep(m, 2, h.k, FUN='/')
dfDat = data.frame(m2[1:90,])

## for each gene fill the missing values
# step 1 = find missing location
# step 2 = simulate new value and return new vector

step1 = function(v) return(which(is.na(v)))

step2 = function(v){
  i = step1(v)
  if (length(i) < 1 || length(i) == length(v)) return(v)
  v.new = replicate(length(i), simpost(na.omit(v)))
  v[i] = v.new
  return(v)
}

dfDat.new = data.frame()

for(i in 1:nrow(dfDat)){
  g1 = as.numeric(dfDat[i,])
  dat.new = tapply(g1, fGroups, step2)
  dfDat.new = rbind(dfDat.new, c(dat.new[[2]], dat.new[[1]]))
}

dimnames(dfDat.new) = dimnames(dfDat)

# # draw a heatmap of the data
# library(NMF)
# source('../CGraphClust/CGraphClust.R')
# 
# m1 = as.matrix(dfDat.new)
# colnames(m1) = fGroups
# m1 = na.omit(m1)
# 
# mCounts = t(m1)
# mCounts = scale(mCounts)
# mCounts = t(mCounts)
# # threshhold the values
# mCounts[mCounts < -3.5] = -3.5
# mCounts[mCounts > 3.5] = 3.5
# 
# # draw the heatmap  color='-RdBu:50'
# aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
#          annColors=NA, Colv=NA)
# 
# # repeat after smoothing
# m1 = t(apply(m1, 1, function(x) f_ivStabilizeData(x, fGroups)))
# colnames(m1) = fGroups
# 
# mCounts = t(m1)
# mCounts = scale(mCounts)
# mCounts = t(mCounts)
# # threshhold the values
# mCounts[mCounts < -3.5] = -3.5
# mCounts[mCounts > 3.5] = 3.5
# 
# # draw the heatmap  color='-RdBu:50'
# aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
#          annColors=NA, Colv=NA)


m1 = as.matrix(dfDat.new)
colnames(m1) = fGroups
m1 = na.omit(m1)

p = apply(m1, 1, function(x) t.test(x ~ fGroups)$p.value)
i = which(p < 0.05)

mCounts = (m1[i,])
mCounts = t(scale(mCounts))

pr.out = prcomp(mCounts, scale=F)
biplot(pr.out)

d = dist(mCounts)
hc = hclust(d)
plot(hc)




