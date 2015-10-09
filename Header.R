# Name: Header.R
# Auth: u.niazi@imperial.ac.uk
# Date: 26 Aug 2015
# Desc: header file with functions and libraries to load

######################### libraries 
library(annotate)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(lumi)
library(limma)
library(downloader)
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
library(NMF)
source('../CGraphClust/CGraphClust.R')
######################### end libraries


######################## global variables
p.old = par()

########################


######################### functions
# Name: f_Plot3DPCA
# Args: mComp = n X 3 matrix with first 3 components as the 3 column vectors
#       color = colours for the points
#       ... additional arguments to the plot function
# Rets: none
# Desc: takes the first 3 components and plots the data in a 3d plot
f_Plot3DPCA = function(mComp, color, ...) {
  x = mComp[,1]
  y = mComp[,2]
  z = mComp[,3]
  if (!require(scatterplot3d)) stop('scatterplot3d library required')
  scatterplot3d(x, y, z, color, ...)
}


# Function: f_lGetPCAClusterCount
# Desc: Takes first two components of the PCA and counts possible clusters. The function does this by 
#       binning the vector of data into X bins, assigning a class label to each bin, counting how many
#       observations in each bin, any total number of bins with at least one observations. This is calculated
#       for both the components of the pca matrix, and the max number of bins with at least one observation, in
#       first or second dimension is reported along with a data.frame with cluster labels.
# Args: pr.out = principal component object returned by prcomp function
# Rets: returns list with 2 elements: 
#       1 - cluster.count = possible number of clusters in the data
#       2 - cluster.label = data.frame with cluster labels
f_lGetPCAClusterCount = function(pr.out){
  # how many clusters in data, using first 2 components
  x1 = pr.out$x[,1]
  x2 = pr.out$x[,2]
  # bin the data from the 2 components
  h1 = hist(x1, plot=F)
  # give a class label to each bin
  c1 = cut(x1, h1$breaks, labels = 1:(length(h1$mids)))
  h2 = hist(x2, plot=F)
  c2 = cut(x2, h2$breaks, labels = 1:(length(h2$mids)))
  # labels for vectors and the class labels
  dfClust = data.frame(lab=names(x1), c1, c2)
  # get contingency table
  mClust = as.matrix(table(c1 = dfClust$c1, c2 = dfClust$c2))
  # count the max of row and col sums that are not zero
  ir = length(which(rowSums(mClust) != 0))
  ic = length(which(colSums(mClust) != 0))
  iClust.count = ifelse(ir > ic, ir, ic)
  lRet = list(cluster.count=iClust.count, cluster.label=dfClust)
  return(lRet)
}



######################### end functions