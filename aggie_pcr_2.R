# Name: aggie_pcr.R
# Date: 15/12/15
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Desc: analysis of pcr data


# data input
dfDat = read.csv('Data_external/aggie_pcr/input2.csv', header=T, stringsAsFactors = F)

fStatus = colnames(dfDat)
fStatus = fStatus[-1]
fStatus = gsub('(\\w+).\\d+', replacement = '\\1', fStatus)
fStatus = factor(fStatus, levels = c('Noninfected', 'Infected'))
table(fStatus)
head(fStatus)

fDays = dfDat[1,-1]
fDays = gsub('\\w+ \\w(\\d+)$', '\\1', fDays[1,])
fDays = factor(fDays, levels=c('0', '7', '28'))
table(fDays)
head(fDays)

xtabs(~ fStatus + fDays)
fGroups = factor(paste0(fStatus, fDays))
ftable(fDays, fStatus, fGroups)

# assign unique names to genes
dfDat = dfDat[-1,]
rn = paste(dfDat$X, 1:nrow(dfDat), sep='.')
rownames(dfDat) = rn
dfDat = dfDat[,-1]

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

# divide by housekeeping gene
m = as.matrix(dfDat)
m = t(apply(m, 1, as.numeric))
h.k = m['GAPDH.87',]
m2 = sweep(m, 2, h.k, FUN='/')
dfDat = data.frame(m2)

## for each gene fill the missing values
# step 1 = find missing location
# step 2 = simulate new value and return new vector

step1 = function(v) return(which(is.na(v)))

step2 = function(v){
  i = step1(v)
  if (length(i) < 1 || length(i) >= length(v)-1) return(v)
  v.new = replicate(length(i), simpost(na.omit(v)))
  v[i] = v.new
  return(v)
}

dfDat.new = data.frame(fGroups)

for(i in 1:nrow(dfDat)){
  g1 = as.numeric(dfDat[i,])
  dat.new = rep(NA, length.out = nrow(dfDat.new))
  l = levels(fGroups)
  for (r in seq_along(l)){
    pos = which(fGroups == l[r])
    dat.new[pos] = step2(g1[pos])
  }
  #dat.new = tapply(g1, fGroups, step2)
  dfDat.new = cbind(dfDat.new, dat.new)
}

colnames(dfDat.new) = c('fGroups', rownames(dfDat))



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


# any of the p.values significant under anova
m1 = as.matrix(dfDat.new[,2:ncol(dfDat.new)])
p.ano = apply(m1, 2, function(x) anova(lm(x ~ fGroups))$Pr[1])
p.adj = p.adjust(p.ano, 'bonf')
table((p.adj < 0.05))

## pairwise t.test
p = apply(m1, 2, function(x) {
  df = data.frame(x, fGroups)
  df = na.omit(df)
  # check if any particular group is too small
  fr = which(table(df$fGroups) < 2)
  i = which(df$fGroups %in% names(fr))
  p.val = (pairwise.t.test(df$x[-i], df$fGroups[-i], p.adjust.method = 'bonf')$p.value)
  return(any(p.val < 0.05))
})
table(p)
i = which(p == TRUE)
x = m1[,i]
df = data.frame(x, fGroups)
df = na.omit(df)

boxplot(x ~ fGroups, main=names(i))
table(df$fGroups)

## PCA plot for the data
f = colSums(m1)
f = is.na(f)
# remove data with NA
mCounts = m1[,!f]
# scale the data
mCounts = scale(mCounts)
f = colSums(mCounts)
f = is.finite(f)
# remove data with NAN
mCounts = mCounts[,f]




m = na.omit(scale(mCounts))
pr.out = prcomp(m, scale=F)
biplot(pr.out)
mCounts = t(scale(mCounts))

pr.out = prcomp(mCounts, scale=F)
biplot(pr.out)

d = dist(mCounts)
hc = hclust(d)
plot(hc)




